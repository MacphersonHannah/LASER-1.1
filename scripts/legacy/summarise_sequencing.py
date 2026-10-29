#!/usr/bin/env python3
"""
Summarise long-read RNA-seq outputs for Results intro figure.

Outputs:
  - per_sample_reads.tsv
  - isoform_summary.tsv
  - (optional) junction_summary.tsv if SQANTI has 'junction_category'
  - summary.json  (rolled-up metrics)

Inputs:
  - --samples CSV with columns:
        sample,bam[,fastq,library]
    (library is optional; falls back to sample if missing)
  - --sqanti  path to SQANTI classification file (tsv/txt)
  - --roi     genomic region, eg 'chr4:90640000-90710000'  (optional)
  - --gene    gene symbol, eg 'SNCA'                       (optional)
  - --min-mapq INTEGER (default 0)
  - --primary-only (default true) count only primary alignments
  - --use-fastq-total  count raw reads from FASTQ for total_reads (slow on very large files)
  - --outdir  output directory

Notes:
  - Requires: pandas, pysam (and gzip if FASTQs are .gz)
  - SQANTI classification is assumed to have columns like:
      'isoform', 'structural_category', and a gene column among:
      ['gene', 'gene_name', 'associated_gene', 'Gene', 'GeneName'].
"""

import argparse, os, sys, json, gzip, re
from pathlib import Path
import pandas as pd
import pysam

GENE_COL_CANDIDATES = ["gene", "gene_name", "associated_gene", "Gene", "GeneName"]

def parse_args():
    ap = argparse.ArgumentParser(description="Summarise long-read RNA-seq outputs.")
    ap.add_argument("--samples", required=True, help="CSV with columns: sample,bam[,fastq,library]")
    ap.add_argument("--sqanti", required=True, help="SQANTI classification file (tsv/txt)")
    ap.add_argument("--roi", default=None, help="Region of interest, eg chr4:90640000-90710000")
    ap.add_argument("--gene", default=None, help="Gene symbol of interest, eg SNCA")
    ap.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ to count as mapped")
    ap.add_argument("--primary-only", action="store_true", default=True, help="Count only primary alignments (default)")
    ap.add_argument("--no-primary-only", dest="primary_only", action="store_false", help="Include secondary/supplementary")
    ap.add_argument("--use-fastq-total", action="store_true", help="Compute total_reads from FASTQ (slow)")
    ap.add_argument("--outdir", default="seq_summary", help="Output directory")
    return ap.parse_args()

def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)

def read_samples_csv(path_csv):
    df = pd.read_csv(path_csv)
    req = {"sample","bam"}
    missing = req - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] --samples is missing required columns: {missing}")
    if "library" not in df.columns:
        df["library"] = df["sample"]
    if "fastq" not in df.columns:
        df["fastq"] = pd.NA
    return df

def count_fastq_reads(fq_path):
    if pd.isna(fq_path) or fq_path == "" or not Path(fq_path).exists():
        return None
    opener = gzip.open if str(fq_path).endswith(".gz") else open
    n = 0
    with opener(fq_path, "rt") as fh:
        for _ in fh:
            n += 1
    return n // 4  # 4 lines per read

def n50(lengths):
    if not lengths:
        return None
    arr = sorted(lengths, reverse=True)
    half = sum(arr) / 2
    run = 0
    for L in arr:
        run += L
        if run >= half:
            return L
    return arr[-1]

def parse_roi(roi_str):
    if not roi_str:
        return None
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", roi_str.replace(",", ""))
    if not m:
        sys.exit(f"[ERROR] Bad --roi '{roi_str}'. Use chr:start-end")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if end <= start:
        sys.exit(f"[ERROR] ROI end must be > start: {roi_str}")
    return chrom, start, end

def count_bam_metrics(bam_path, min_mapq=0, primary_only=True, roi=None):
    """
    Returns dict with:
      total_primary_aligns, mapped_primary_aligns, roi_primary_aligns,
      mean_read_len, n50_read_len, median_mapq
    NB 'total_primary_aligns' counts primary alignments present in the BAM.
    """
    if not Path(bam_path).exists():
        sys.exit(f"[ERROR] BAM not found: {bam_path}")
    bam = pysam.AlignmentFile(bam_path, "rb")

    # ensure index for ROI
    if roi and not Path(bam_path + ".bai").exists() and not Path(bam_path.replace(".bam",".bai")).exists():
        sys.exit(f"[ERROR] BAM index (.bai) required for ROI fetch: {bam_path}")

    def is_primary(aln):
        return not (aln.is_secondary or aln.is_supplementary)

    total = 0
    mapped = 0
    mapqs = []
    read_lengths = []

    # iterate entire BAM
    for aln in bam.fetch(until_eof=True):
        if primary_only and not is_primary(aln):
            continue
        total += 1
        if not aln.is_unmapped and aln.mapping_quality >= min_mapq:
            mapped += 1
            mapqs.append(aln.mapping_quality)
        # read length (query length; for long-read RNA this is useful)
        ql = aln.query_length if aln.query_length is not None else 0
        if ql:
            read_lengths.append(ql)

    roi_count = None
    if roi:
        chrom, start, end = roi
        roi_count = 0
        for aln in bam.fetch(chrom, start, end):
            if primary_only and not is_primary(aln):
                continue
            if aln.is_unmapped or aln.mapping_quality < min_mapq:
                continue
            # any overlap counts
            roi_count += 1

    bam.close()
    # summarise
    mean_len = (sum(read_lengths)/len(read_lengths)) if read_lengths else None
    n50_len = n50(read_lengths)
    median_mapq = (int(pd.Series(mapqs).median()) if mapqs else None)
    return dict(
        total_primary_aligns=total,
        mapped_primary_aligns=mapped,
        roi_primary_aligns=roi_count,
        mean_read_len=mean_len,
        n50_read_len=n50_len,
        median_mapq=median_mapq,
    )

def load_sqanti(path_sq):
    # SQANTI files are usually tab-separated
    df = pd.read_csv(path_sq, sep=None, engine="python")
    # normalise colnames
    df.columns = [c.strip() for c in df.columns]
    # find gene column
    gene_col = None
    for c in GENE_COL_CANDIDATES:
        if c in df.columns:
            gene_col = c
            break
    if gene_col is None:
        # tolerate absence; gene-based summaries will be NA
        gene_col = None
    if "isoform" not in df.columns:
        # TALON outputs may use 'transcript' or similar
        for cand in ["transcript","pbid","isoform_id","isoformName"]:
            if cand in df.columns:
                df = df.rename(columns={cand:"isoform"})
                break
    if "isoform" not in df.columns:
        sys.exit("[ERROR] Could not find an isoform identifier column in SQANTI file.")
    if "structural_category" not in df.columns:
        # SQANTI3 uses 'structural_category'; older may use 'category'
        if "category" in df.columns:
            df = df.rename(columns={"category":"structural_category"})
        else:
            # will skip category summaries if missing
            df["structural_category"] = pd.NA
    return df, gene_col

def summarise_isoforms(df_sq, gene_col=None, gene=None):
    # total unique isoforms overall
    total_iso = df_sq["isoform"].nunique()

    # gene filter if possible
    if gene and gene_col and gene_col in df_sq.columns:
        df_gene = df_sq[df_sq[gene_col].astype(str) == gene]
        gene_iso = df_gene["isoform"].nunique()
    else:
        df_gene = pd.DataFrame(columns=df_sq.columns)
        gene_iso = None

    # structural category tallies
    cat_all = (
        df_sq.dropna(subset=["structural_category"])
             .groupby("structural_category")["isoform"]
             .nunique()
             .reset_index(name="all_isoforms")
    )
    if not df_gene.empty:
        cat_gene = (
            df_gene.dropna(subset=["structural_category"])
                  .groupby("structural_category")["isoform"]
                  .nunique()
                  .reset_index(name="gene_isoforms")
        )
        cats = pd.merge(cat_all, cat_gene, on="structural_category", how="outer").fillna(0).sort_values("all_isoforms", ascending=False)
    else:
        cats = cat_all
        cats["gene_isoforms"] = 0

    # optional junction summary if present
    junc_df = None
    if "junction_category" in df_sq.columns:
        j_all = (
            df_sq.dropna(subset=["junction_category"])
                 .groupby("junction_category")["isoform"]
                 .nunique()
                 .reset_index(name="all_isoforms")
        )
        if not df_gene.empty:
            j_gene = (
                df_gene.dropna(subset=["junction_category"])
                      .groupby("junction_category")["isoform"]
                      .nunique()
                      .reset_index(name="gene_isoforms")
            )
            junc_df = pd.merge(j_all, j_gene, on="junction_category", how="outer").fillna(0).sort_values("all_isoforms", ascending=False)
        else:
            j_all["gene_isoforms"] = 0
            junc_df = j_all

    return dict(
        total_isoforms=int(total_iso),
        gene_isoforms=(int(gene_iso) if gene_iso is not None else None),
        cats=cats,
        junc_df=junc_df
    )

def main():
    a = parse_args()
    ensure_dir(a.outdir)

    # ROI parse (once)
    roi = parse_roi(a.roi) if a.roi else None

    # Samples
    samples = read_samples_csv(a.samples)

    # Per-sample read/mapping summaries
    rows = []
    for _, r in samples.iterrows():
        sample = str(r["sample"])
        library = str(r["library"])
        bam = str(r["bam"])
        fastq = r["fastq"] if "fastq" in r and not pd.isna(r["fastq"]) else None

        metrics = count_bam_metrics(
            bam_path=bam,
            min_mapq=a.min_mapq,
            primary_only=a.primary_only,
            roi=roi
        )
        # totals
        if a.use_fastq_total and fastq:
            total_reads = count_fastq_reads(fastq)
        else:
            total_reads = metrics["total_primary_aligns"]  # approximation when FASTQ not provided

        mapped = metrics["mapped_primary_aligns"]
        mapping_rate = (mapped / total_reads) if total_reads and total_reads > 0 else None
        roi_reads = metrics["roi_primary_aligns"]
        roi_frac = (roi_reads / mapped) if (roi_reads is not None and mapped > 0) else None

        rows.append(dict(
            sample=sample,
            library=library,
            bam=bam,
            fastq=(fastq if fastq else ""),
            total_reads=(int(total_reads) if total_reads is not None else None),
            mapped_reads=int(mapped),
            mapping_rate=(round(mapping_rate, 4) if mapping_rate is not None else None),
            roi_reads=(int(roi_reads) if roi_reads is not None else None),
            roi_fraction=(round(roi_frac, 4) if roi_frac is not None else None),
            mean_read_len=(round(metrics["mean_read_len"], 1) if metrics["mean_read_len"] else None),
            n50_read_len=metrics["n50_read_len"],
            median_mapq=metrics["median_mapq"]
        ))

    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(Path(a.outdir) / "per_sample_reads.tsv", sep="\t", index=False)

    # SQANTI summaries
    sq, gene_col = load_sqanti(a.sqanti)
    iso_sum = summarise_isoforms(sq, gene_col=gene_col, gene=a.gene)

    # Isoform category summary table
    cats = iso_sum["cats"].copy()
    cats.to_csv(Path(a.outdir) / "isoform_summary.tsv", sep="\t", index=False)

    # Optional junction summary
    if iso_sum["junc_df"] is not None:
        iso_sum["junc_df"].to_csv(Path(a.outdir) / "junction_summary.tsv", sep="\t", index=False)

    # JSON rollup
    rollup = {
        "n_samples": int(per_sample["sample"].nunique()),
        "roi": a.roi,
        "gene": a.gene,
        "min_mapq": a.min_mapq,
        "primary_only": a.primary_only,
        "totals": {
            "reads_total_sum": int(per_sample["total_reads"].fillna(0).sum()),
            "reads_mapped_sum": int(per_sample["mapped_reads"].sum()),
            "roi_reads_sum": (int(per_sample["roi_reads"].fillna(0).sum()) if "roi_reads" in per_sample.columns else None),
        },
        "isoforms": {
            "total_unique_isoforms": iso_sum["total_isoforms"],
            "gene_unique_isoforms": iso_sum["gene_isoforms"],
            "by_structural_category": cats.to_dict(orient="records"),
        }
    }
    with open(Path(a.outdir) / "summary.json", "w") as fh:
        json.dump(rollup, fh, indent=2)

    # stdout summary
    print("\n=== Sequencing summary ===")
    print(per_sample.head().to_string(index=False))
    print("\nIsoforms (overall):", iso_sum["total_isoforms"])
    if a.gene:
        print(f"Isoforms ({a.gene}):", iso_sum["gene_isoforms"])
    if a.roi:
        print(f"ROI reads total: {int(per_sample['roi_reads'].fillna(0).sum())}")
    print(f"\nTables written to: {a.outdir}\n")

if __name__ == "__main__":
    main()
