#!/usr/bin/env python3

## 14/11/2025 added ability to use indels, not just SNVs
## 05/12/2025 Changed the logic to improve specificity but still include small unassigned (by Sqanti) reads. For each gene (as defined by the input bed)
## 1) take all reads that intersect the genomic region
## 2) remove any that are in the wrong direction
## 3) remove any that Sqanti has assigned to a different gene (except fusions of gene of interest)
## Some reads may remain unassigned by Sqanti (no transcript match); these are kept if they overlap the region and are in the right direction. May need manual checking by user.

"""
Fast version (15/10/2025) - same output, optimised for speed.

Safe optimisations:
  1. Fetch only merged regions from BAM (no full-contig iteration)
  2. Cache variants in memory for DeepVariant + WhatsHap VCFs
  3. Parallelise by chromosome (each process opens its own files)

Now tracks:
  - SNPs vs non-SNP variants (indels + MNPs)
  - total phased variants per region
and uses total phased variants for ASE QC.
"""

import argparse, os, sys, csv
from typing import List, Dict, Any, Tuple
import pysam
from concurrent.futures import ProcessPoolExecutor
import sqlite3
import re

import numpy as np

 
# Helpers

def ensure_index(vcf_or_bam_path: str, ext: str) -> str:
    """Ensure a file (VCF or BAM) is indexed."""
    if ext == "vcf":
        if os.path.exists(vcf_or_bam_path + ".tbi") or os.path.exists(vcf_or_bam_path + ".csi"):
            return vcf_or_bam_path
        if not vcf_or_bam_path.endswith(".vcf.gz"):
            sys.stderr.write(f"[ERROR] Must be bgzipped: {vcf_or_bam_path}\n")
            sys.exit(1)
        pysam.tabix_index(vcf_or_bam_path, preset="vcf", force=True)
    elif ext == "bam":
        if not os.path.exists(vcf_or_bam_path + ".bai"):
            pysam.index(vcf_or_bam_path)
    return vcf_or_bam_path


def load_bed_rows(bed_path: str) -> List[Dict[str, Any]]:
    """
    Expect BED-like file with at least:
      col1: chrom
      col2: start
      col3: end
      col4: gene symbol (optional)
      col5: Ensembl gene ID (optional)
      col6: strand (+/-) (optional)

    If col6 is missing, we assume '+' strand for now. - For testing, remove this.
    """
    regions = []
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            # If strand column is missing, assume '+'
            strand = "+"
            if len(parts) > 5 and parts[5] in ("+", "-"):
                strand = parts[5]

            regions.append({
                "chrom": parts[0],
                "start": int(float(parts[1])),
                "end": int(float(parts[2])),
                "gene": parts[3] if len(parts) > 3 else "",
                "ensembl_id": parts[4] if len(parts) > 4 else "",
                "strand": strand,
            })
    return regions



def merge_intervals(regions: List[Dict[str, Any]], chrom: str) -> List[Tuple[int, int]]:
    """Merge overlapping or adjacent intervals per chromosome."""
    ivals = sorted([(r["start"], r["end"]) for r in regions if r["chrom"] == chrom])
    merged = []
    for s, e in ivals:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return merged


def fmt_props(h1, h2, un):
    tot = h1 + h2 + un
    if tot == 0:
        return "0:0:0"
    p1 = int(round(100.0 * h1 / tot))
    p2 = int(round(100.0 * h2 / tot))
    pu = max(0, 100 - p1 - p2)
    return f"{p1}:{p2}:{pu}"


def is_snp(rec) -> bool:
    """True if variant is a 1bp ref and 1bp alt (SNP)."""
    return bool(rec.alts) and len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts)


def is_non_snp_variant(rec) -> bool:
    """
    Non-SNP variant: anything that is not a pure 1bp ref/alt SNP.
    This includes indels and MNPs. 
    """
    return not is_snp(rec)


def build_bin_index(regions, bin_size=100_000):
    idx = {}
    for i, r in enumerate(regions):
        chrom, start, end = r["chrom"], r["start"], r["end"]
        if end <= start:
            continue
        d = idx.setdefault(chrom, {})
        b0 = start // bin_size
        b1 = (end - 1) // bin_size
        for b in range(b0, b1 + 1):
            d.setdefault(b, []).append(i)
    return idx


def overlapping_region_indices(bin_index, chrom, qstart, qend, bin_size=100_000):
    d = bin_index.get(chrom)
    if not d:
        return []
    b0 = qstart // bin_size
    b1 = max(b0, (max(qstart, qend - 1)) // bin_size)
    out = []
    for b in range(b0, b1 + 1):
        if b in d:
            out.extend(d[b])
    return out

 
# Pre-cache VCFs (variants)
 

def load_variants_dv(vcf_path, sample_name, only_pass, min_minor_count, min_minor_frac):
    """
    Cache positions of usable variants per chrom for DeepVariant VCF.
    Splits into SNPs vs non-SNPs (indels + MNPs).
    """
    vcf = pysam.VariantFile(vcf_path)
    snps = {}
    non_snps = {}

    for chrom in vcf.header.contigs:
        snps[chrom] = []
        non_snps[chrom] = []

        for rec in vcf.fetch(chrom):
            if only_pass and rec.filter and "PASS" not in rec.filter.keys():
                continue

            # Optional AD-based filter (for DV sample)
            if sample_name and sample_name in rec.samples:
                call = rec.samples[sample_name]
                AD = call.get("AD", None)
                if AD and len(AD) >= 2:
                    ref = int(AD[0])
                    alt = int(max(AD[1:]))
                    tot = ref + alt
                    if tot > 0:
                        minor = min(ref, alt)
                        if minor < min_minor_count or (minor / tot) < min_minor_frac:
                            continue

            if is_snp(rec):
                snps[chrom].append(rec.pos)
            elif is_non_snp_variant(rec):
                non_snps[chrom].append(rec.pos)

    vcf.close()
    return snps, non_snps


def load_variants_wh(vcf_path, sample_name):
    """
    Cache positions of phased heterozygous variants per chrom for WhatsHap VCF.
    Splits into SNPs vs non-SNPs (indels + MNPs).
    """
    vcf = pysam.VariantFile(vcf_path)
    snps = {}
    non_snps = {}

    for chrom in vcf.header.contigs:
        snps[chrom] = []
        non_snps[chrom] = []

        for rec in vcf.fetch(chrom):
            call = vcf.header.samples and rec.samples.get(sample_name)
            if not call:
                continue

            gt = call.get("GT", None)
            is_phased = getattr(call, "phased", False)

            if not gt or len(gt) != 2 or gt[0] == gt[1] or not is_phased:
                continue

            if is_snp(rec):
                snps[chrom].append(rec.pos)
            elif is_non_snp_variant(rec):
                non_snps[chrom].append(rec.pos)

    vcf.close()
    return snps, non_snps


def load_variants_wh_full(vcf_path, sample_name):
    """
    Cache phased heterozygous WhatsHap variants per chrom with full info:
    pos, ref, alt (comma-joined), id, is_snp (1/0).
    only for optional export
    """
    vcf = pysam.VariantFile(vcf_path)
    tmp = {}

    for chrom in vcf.header.contigs:
        tmp[chrom] = []
        for rec in vcf.fetch(chrom):
            call = vcf.header.samples and rec.samples.get(sample_name)
            if not call:
                continue

            gt = call.get("GT", None)
            is_phased = getattr(call, "phased", False)

            # EXACTLY the same variant inclusion logic as load_variants_wh:
            # phased + heterozygous only
            if not gt or len(gt) != 2 or gt[0] == gt[1] or not is_phased:
                continue

            alt = ",".join(rec.alts) if rec.alts else ""
            vid = rec.id if rec.id and rec.id != "." else ""
                        # --- Variant QC fields (from VCF record + sample call) ---
            qual = rec.qual if rec.qual is not None else ""
            filt = ";".join(list(rec.filter.keys())) if rec.filter is not None else ""
            if filt == "":
                filt = "PASS"

            # FORMAT fields from the sample call
            gt_str = ""
            ps = ""
            dp = ""
            gq = ""
            vaf = ""

            gt_val = call.get("GT", None)
            if gt_val is not None and len(gt_val) == 2:
                sep = "|" if is_phased else "/"
                gt_str = f"{gt_val[0]}{sep}{gt_val[1]}"

            ps_val = call.get("PS", None)
            if ps_val is not None:
                ps = str(ps_val)

            dp_val = call.get("DP", None)
            if dp_val is not None:
                dp = str(dp_val)

            gq_val = call.get("GQ", None)
            if gq_val is not None:
                gq = str(gq_val)

            vaf_val = call.get("VAF", None)
            # VAF is Number=A so may be list/tuple
            if vaf_val is not None:
                if isinstance(vaf_val, (list, tuple)):
                    vaf = ",".join(str(x) for x in vaf_val)
                else:
                    vaf = str(vaf_val)

            tmp[chrom].append((
                int(rec.pos),
                rec.ref,
                alt,
                vid,
                1 if is_snp(rec) else 0,
                filt,
                qual,
                gt_str,
                ps,
                dp,
                gq,
                vaf
            ))


    vcf.close()

    # Pack into sorted arrays for fast slicing
    packed = {}
    for chrom, lst in tmp.items():
        # EMPTY chromosome: store empty arrays/lists
        if not lst:
            packed[chrom] = {
                "pos": np.array([], dtype=np.int32),
                "ref": [],
                "alt": [],
                "id": [],
                "is_snp": np.array([], dtype=np.int8),
                "filter": [],
                "qual": [],
                "gt": [],
                "ps": [],
                "dp": [],
                "gq": [],
                "vaf": [],
            }
            continue

        # NON-empty: sort + pack real values
        lst.sort(key=lambda x: x[0])
        packed[chrom] = {
            "pos": np.array([x[0] for x in lst], dtype=np.int32),
            "ref": [x[1] for x in lst],
            "alt": [x[2] for x in lst],
            "id":  [x[3] for x in lst],
            "is_snp": np.array([x[4] for x in lst], dtype=np.int8),
            "filter": [x[5] for x in lst],
            "qual": [x[6] for x in lst],
            "gt": [x[7] for x in lst],
            "ps": [x[8] for x in lst],
            "dp": [x[9] for x in lst],
            "gq": [x[10] for x in lst],
            "vaf": [x[11] for x in lst],
        }

    return packed



def slice_variants_wh_full(varpack, chrom, start, end):
    """
    Return list of (pos, ref, alt, id, is_snp) for variants in [start, end).
    """
    d = varpack.get(chrom)
    if not d:
        return []
    pos = d["pos"]
    if pos.size == 0:
        return []

    left = np.searchsorted(pos, start, side="left")
    right = np.searchsorted(pos, end, side="left")

    out = []
    for i in range(left, right):
        out.append((
            int(pos[i]),
            d["ref"][i],
            d["alt"][i],
            d["id"][i],
            int(d["is_snp"][i]),
            d["filter"][i],
            d["qual"][i],
            d["gt"][i],
            d["ps"][i],
            d["dp"][i],
            d["gq"][i],
            d["vaf"][i],
        ))
    return out



def to_numpy_positions(pos_dict):
    """Convert per-chrom lists of positions to NumPy arrays for fast searchsorted."""
    for chrom in list(pos_dict.keys()):
        pos_dict[chrom] = np.array(pos_dict[chrom], dtype=np.int32)
    return pos_dict


def count_variants_cached(pos_dict, chrom, start, end):
    arr = pos_dict.get(chrom)
    if arr is None or arr.size == 0:
        return 0
    left = np.searchsorted(arr, start, side="left")
    right = np.searchsorted(arr, end, side="left")
    return int(right - left)

 
# Worker per chromosome
 

def process_chrom(args_tuple):
    (
        chrom,
        chrom_regions,
        dv_snps,
        dv_non_snps,
        ph_snps,
        ph_non_snps,
        bam_path,
        args,
        read_to_gene,
        ph_varpack,
    ) = args_tuple

    rows = []
    var_rows = []
    merged_intervals = merge_intervals(chrom_regions, chrom)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # pre-count haplotagged reads in all regions
        h1 = [0] * len(chrom_regions)
        h2 = [0] * len(chrom_regions)
        un = [0] * len(chrom_regions)

        # build per-chrom bin index once
        bin_index = build_bin_index(chrom_regions, bin_size=100_000)

        for s, e in merged_intervals:
            for aln in bam.fetch(chrom, s, e):
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue

                rstart = aln.reference_start
                rend = aln.reference_end or (rstart + aln.query_length)
                tag = aln.get_tag(args.hp_tag) if aln.has_tag(args.hp_tag) else None

                # Read strand from BAM flags
                read_strand = '-' if aln.is_reverse else '+'

                # SQANTI-associated gene string for this read (can be "", Ensembl, fusions, novelGene_xxx, etc.)
                gene_str = read_to_gene.get(aln.query_name, "")  # default to empty string

                cand = overlapping_region_indices(bin_index, chrom, rstart, rend, bin_size=100_000)
                if not cand:
                    continue

                for idx in cand:
                    R = chrom_regions[idx]
                    if rend <= R["start"] or rstart >= R["end"]:
                        continue

                    # 1) Directionality: require read strand matches gene strand.
                    #    For now, gene strand defaults to '+' if BED didn't provide it.
                    gene_strand = R.get("strand", "+") or "+"
                    if gene_strand in ('+', '-') and read_strand != gene_strand:
                        continue

                    # 2) SQANTI gene filter:
                    #    - If gene_str is non-empty: enforce Ensembl match (including fusions)
                    #    - If gene_str is empty:
                    #         * include only if args.include_unassigned is True

                    # BED: Ensembl gene ID (likely without version)
                    bed_gene_ens = (R.get("ensembl_id") or "").strip()
                    if not bed_gene_ens:
                        # If the BED doesn't have an Ensembl ID for this region, be strict and skip.
                        continue

                    bed_gene_ens_base = bed_gene_ens.split('.')[0]

                    if gene_str:
                        # SQANTI-associated gene string: Ensembl ID or fusion of IDs
                        # eg "ENSG00000145335.17" or "ENSG00000145335.17_ENSG000002XXXX"
                        tokens = re.split(r'[_;,|\s]+', gene_str)
                        tokens_clean = [t.split('.')[0] for t in tokens if t]
                        ok = bed_gene_ens_base in tokens_clean


                        ok = bed_gene_ens_base in tokens_clean
                        if not ok:
                            # SQANTI says this read belongs to some other gene(s); skip for this BED locus
                            continue
                        # else: read is assigned to this gene (or a fusion including it)
                    else:
                        # gene_str == "" -> no SQANTI assignment.
                        # Only keep if the user explicitly asked to include unassigned reads.
                        if not args.include_unassigned:
                            continue
                        # if include_unassigned: fall through and keep this read


                    # If we reach here, the read is considered to belong to THIS locus
                    if tag == 1:
                        h1[idx] += 1
                    elif tag == 2:
                        h2[idx] += 1
                    else:
                        un[idx] += 1

    # Now QC per region
    sample_label = os.path.basename(args.out).replace(".loci.csv", "").replace(".csv", "")

    for i, R in enumerate(chrom_regions):
        total_assigned = h1[i] + h2[i]
        if total_assigned < args.min_total_reads:
            status = "TOO_FEW_TOTAL_READS"
        else:
            status = None  # continue to variant QC

        # DV variant counts
        dv_snps_count = count_variants_cached(dv_snps, chrom, R["start"], R["end"])
        dv_non_snps_count = count_variants_cached(dv_non_snps, chrom, R["start"], R["end"])

        # Phased variant counts (WhatsHap)
        phased_snps = count_variants_cached(ph_snps, chrom, R["start"], R["end"])
        phased_non_snps = count_variants_cached(ph_non_snps, chrom, R["start"], R["end"])
        total_phased_variants = phased_snps + phased_non_snps

        variant_qc = "PASS" if total_phased_variants >= args.min_variant_sites else "FAIL"

        if variant_qc == "FAIL":
            status = "NOT_ENOUGH_PHASED_VARIANTS"
        else:
            tot = h1[i] + h2[i] + un[i]
            low = (h1[i] < args.min_reads_per_haplotype) or (h2[i] < args.min_reads_per_haplotype)
            un_ok = not (tot > 0 and (un[i] / tot) > args.max_unphased_frac)
            if not low and un_ok:
                status = "PASS"
            elif low and not un_ok:
                status = "FEW_HAP_READS_AND_HIGH_UNPHASED"
            elif low:
                status = "TOO_FEW_HAPLOTYPE_READS"
            else:
                status = "TOO_MANY_UNPHASED_READS"

        if ph_varpack is not None:
            vars_in_locus = slice_variants_wh_full(ph_varpack, chrom, R["start"], R["end"])
            for (vpos, vref, valt, vid, vis_snp, vfilt, vqual, vgt, vps, vdp, vgq, vvaf) in vars_in_locus:
                var_rows.append([
                    sample_label,
                    chrom,
                    R["start"],
                    R["end"],
                    R["gene"],
                    R["ensembl_id"],
                    variant_qc,  # locus-level "enough variants" flag; keep if you want
                    vpos,
                    vref,
                    valt,
                    vid,
                    "SNP" if vis_snp == 1 else "nonSNP",
                    vfilt,
                    vqual,
                    vgt,
                    vps,
                    vdp,
                    vgq,
                    vvaf,
                ])



        rows.append([
            sample_label,
            chrom,
            R["start"],
            R["end"],
            R["gene"],
            R["ensembl_id"],
            dv_snps_count,
            dv_non_snps_count,
            phased_snps,
            phased_non_snps,
            total_phased_variants,
            variant_qc,
            h1[i],
            h2[i],
            un[i],
            status,
        ])

    return rows, var_rows




## Adding the mapping from the SQlite (per read, which gene does SQANTI assign it to).Instead of only using genomic coordinates.
## this one works via SQanti categorisation
# def load_read_gene_map(sqlite_path: str, sample_label: str) -> Dict[str, str]:
#     """
#     Build mapping: read_id -> sqanti_associated_gene string (Ensembl gene ID(s))
#     for a given sample.

#     - joins reads.pbid -> transcripts.pbid
#     - uses transcripts.sqanti_associated_gene (Ensembl ID; may contain fusions like ENSG1_ENSG2)
#     """
#     conn = sqlite3.connect(sqlite_path)
#     cur = conn.cursor()

#     query = """
#         SELECT r.read_id, t.sqanti_associated_gene
#         FROM reads r
#         JOIN transcripts t ON r.pbid = t.pbid
#         WHERE r.sample = ?
#     """

#     mapping: Dict[str, str] = {}
#     for read_id, gene_str in cur.execute(query, (sample_label,)):
#         if read_id is None:
#             continue
#         mapping[read_id] = gene_str or ""

#     conn.close()
#     return mapping


## this one works via featurecounts characterisation
def load_read_gene_map(sqlite_path: str, sample_label: str) -> Dict[str, str]:
    """
    read_id -> featureCounts gene string for this sample.

    Pull from reads table so sample labels are correct.
    Returns "" if status != Assigned or genes is NULL/NA.
    """
    conn = sqlite3.connect(sqlite_path)
    cur = conn.cursor()

    # Filter to Assigned reads only at the SQL level -- avoids loading millions of
    # unassigned rows into the Python dict. Unassigned reads map to "" anyway.
    query = """
        SELECT read_id, featurecounts_gene
        FROM reads
        WHERE sample = ?
          AND featurecounts_status = 'Assigned'
          AND featurecounts_gene IS NOT NULL
          AND featurecounts_gene != 'NA'
          AND featurecounts_gene != ''
    """

    mapping: Dict[str, str] = {}
    for read_id, genes in cur.execute(query, (sample_label,)):
        if not read_id:
            continue
        mapping[read_id] = genes or ""  # already filtered to Assigned + non-null above

    conn.close()
    return mapping




 
# Main
 

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dv_vcf", required=True)
    ap.add_argument("--phased_vcf", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only_pass", action="store_true")
    # now use total variants (SNPs + indels/MNPs)
    ap.add_argument("--min_variant_sites", type=int, default=1)
    ap.add_argument("--min_total_reads", type=int, default=0)
    ap.add_argument("--min_minor_count", type=int, default=0)
    ap.add_argument("--min_minor_frac", type=float, default=0.0)
    ap.add_argument("--dv_sample", default=None)
    ap.add_argument("--phased_sample", default=None)
    ap.add_argument("--min_reads_per_haplotype", type=int, default=5)
    ap.add_argument("--max_unphased_frac", type=float, default=0.3)
    ap.add_argument("--hp_tag", default="HP")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--sqlite", required=True, help="SQLite DB with read->gene mapping")
    ap.add_argument(
    "--include_unassigned",
    action="store_true",
    help="Include reads with no SQANTI-associated gene in locus counts")
    ap.add_argument(
    "--variants_out",
    default=None,
    help="Optional TSV: per-locus phased het WhatsHap variants (pos/ref/alt) overlapping each locus"
)


    args = ap.parse_args()

    dv_vcf_path = ensure_index(args.dv_vcf, "vcf")
    ph_vcf_path = ensure_index(args.phased_vcf, "vcf")
    bam_path = ensure_index(args.bam, "bam")

    regions = load_bed_rows(args.bed)
    chroms = sorted(set(r["chrom"] for r in regions))

    # sample label from output filename; must match 'reads.sample' in SQLite
    sample_label = os.path.basename(args.out).replace(".loci.csv", "").replace(".csv", "")

    # per-read gene annotation from SQLite (for this sample only)
    read_to_gene = load_read_gene_map(args.sqlite, sample_label)
    sys.stderr.write(f"[INFO] Loaded {len(read_to_gene)} read->gene mappings from {args.sqlite}\n")


    # Pre-cache variants once globally
    with pysam.VariantFile(dv_vcf_path) as dv:
        dv_sample = args.dv_sample or (dv.header.samples[0] if dv.header.samples else None)
    with pysam.VariantFile(ph_vcf_path) as pv:
        ph_sample = args.phased_sample or (pv.header.samples[0] if pv.header.samples else None)

    sys.stderr.write("[INFO] Caching variants...\n")
    dv_snps, dv_non_snps = load_variants_dv(
        dv_vcf_path,
        dv_sample,
        args.only_pass,
        args.min_minor_count,
        args.min_minor_frac,
    )
    ph_snps, ph_non_snps = load_variants_wh(ph_vcf_path, ph_sample)
    ph_varpack = None
    if args.variants_out:
        sys.stderr.write("[INFO] Caching full phased variants for export...\n")
        ph_varpack = load_variants_wh_full(ph_vcf_path, ph_sample)

    dv_snps = to_numpy_positions(dv_snps)
    dv_non_snps = to_numpy_positions(dv_non_snps)
    ph_snps = to_numpy_positions(ph_snps)
    ph_non_snps = to_numpy_positions(ph_non_snps)

    # Parallel per chromosome
    sys.stderr.write("[INFO] Processing chromosomes...\n")
    jobs = [
        (
            chrom,
            [r for r in regions if r["chrom"] == chrom],
            dv_snps,
            dv_non_snps,
            ph_snps,
            ph_non_snps,
            bam_path,
            args,
            read_to_gene,
            ph_varpack,   # <-- ADD THIS
        )
        for chrom in chroms
    ]

    all_var_rows = []

    # Main per-locus output
    with open(args.out, "w", newline="") as outfh:
        writer = csv.writer(outfh)
        writer.writerow([
            "sample",
            "chrom",
            "start",
            "end",
            "gene",
            "ensembl_id",
            "DV_SNPs",
            "DV_nonSNP_variants",
            "phased_SNPs",
            "phased_nonSNP_variants",
            "total_phased_variants",
            "variant_QC",
            "h1_reads",
            "h2_reads",
            "unphased_reads",
            "status",
        ])


        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            for loci_rows, var_rows in ex.map(process_chrom, jobs):
                writer.writerows(loci_rows)
                if args.variants_out:
                    all_var_rows.extend(var_rows)

    if args.variants_out:
        with open(args.variants_out, "w", newline="") as vfh:
            vw = csv.writer(vfh, delimiter="\t")
            vw.writerow([
                "sample", "chrom", "start", "end", "gene", "ensembl_id",
                "variant_QC",
                "variant_pos", "ref", "alt", "id", "variant_class",
                "FILTER", "QUAL", "GT", "PS", "DP", "GQ", "VAF"
            ])
            vw.writerows(all_var_rows)

    # Optional per-chromosome summary
    summary_path = args.out.replace(".loci.csv", ".variant_summary.csv")
    with open(summary_path, "w", newline="") as sfh:
        sw = csv.writer(sfh)
        sw.writerow([
            "chrom",
            "DV_SNPs",
            "DV_nonSNP_variants",
            "phased_SNPs",
            "phased_nonSNP_variants",
            "total_phased_variants",
        ])
        for chrom in chroms:
            dv_s = len(dv_snps.get(chrom, []))
            dv_n = len(dv_non_snps.get(chrom, []))
            ph_s = len(ph_snps.get(chrom, []))
            ph_n = len(ph_non_snps.get(chrom, []))
            sw.writerow([chrom, dv_s, dv_n, ph_s, ph_n, ph_s + ph_n])

    sys.stderr.write("[INFO] Done.\n")


if __name__ == "__main__":
    main()
