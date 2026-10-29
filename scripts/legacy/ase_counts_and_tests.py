#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
from collections import defaultdict
from scipy.stats import binomtest

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tagged-bam", dest="bam", required=True)
    ap.add_argument("--passes", required=True)  # passed_successful_loci.reason.wide.csv
    ap.add_argument("--bed", required=True)     # 4 col BED: chrom start end gene
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()

def read_pass_genes(csv_path, sample, out_path):
    import pandas as pd
    import sys

    df = pd.read_csv(csv_path)

    if 'gene' not in df.columns:
        print(f"[ERROR] {csv_path} is missing 'gene' column", file=sys.stderr)
        sys.exit(1)

    if sample not in df.columns:
        print(f"[INFO] No ASE available for sample '{sample}' -- skipping.", file=sys.stderr)

        # Touch empty output file so Nextflow doesn't crash
        open(out_path, 'w').close()

        sys.exit(0)

    return set(df.loc[df[sample].astype(str).str.upper() == 'PASS', 'gene'])




def read_bed_grouped(bed_path):
    """
    Returns dict: gene -> {chrom -> list of (start,end)}.
    If multiple lines per gene, intervals are later merged per chrom.
    """
    gene_map = defaultdict(lambda: defaultdict(list))
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 4:
                continue
            chrom, start, end, gene = toks[:4]
            gene_map[gene][chrom].append((int(start), int(end)))
    # merge overlapping per chrom
    for gene, by_chr in gene_map.items():
        for chrom, ivals in by_chr.items():
            by_chr[chrom] = merge_intervals(sorted(ivals))
    return gene_map

def merge_intervals(intervals):
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        last = merged[-1]
        if s <= last[1]:  # overlap or touch
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]

def count_in_gene(bam, chrom, intervals):
    """
    Count reads overlapping any of the merged intervals on a chrom. deduplicate by read.query_name so a read spanning two merged
    intervals on the same chrom is counted once.
    """
    seen = set()
    counts = defaultdict(int)

    for start, end in intervals:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rid = read.query_name
            # per gene per chrom dedup
            key = (rid, chrom)
            if key in seen:
                continue
            seen.add(key)

            hp = "unphased"
            if read.has_tag("HP"):
                try:
                    hp = str(read.get_tag("HP"))
                except Exception:
                    hp = "unphased"

            # determine if this read was annotated at all (by our tagger)
            has_any_anno = any(read.has_tag(t) for t in ("IN","CD","NM","OR"))
            if not has_any_anno:
                counts[(hp, "UNANNOTATED")] += 1
                continue

            # pull each tag with a fallback
            cd = read.get_tag("CD") if read.has_tag("CD") else "read_unannotated"
            inn = read.get_tag("IN") if read.has_tag("IN") else "read_unannotated"
            nm = read.get_tag("NM") if read.has_tag("NM") else "read_unannotated"
            orf = read.get_tag("OR") if read.has_tag("OR") else "read_unannotated"

            counts[(hp, f"CD:{cd}")]  += 1
            counts[(hp, f"IN:{inn}")] += 1
            counts[(hp, f"NM:{nm}")]  += 1
            counts[(hp, f"OR:{orf}")] += 1
            counts[(hp, "TOTAL")]     += 1
    return counts

def summarise_counts(sample, gene, counts):
    """
    Flatten to a single row with stable column names and a binomial test.
    """
    # totals
    h1 = counts.get(("1","TOTAL"), 0)
    h2 = counts.get(("2","TOTAL"), 0)
    unph = counts.get(("unphased","TOTAL"), 0)
    unanno_h1 = counts.get(("1","UNANNOTATED"), 0)
    unanno_h2 = counts.get(("2","UNANNOTATED"), 0)
    unanno_unph = counts.get(("unphased","UNANNOTATED"), 0)

    total_phased = h1 + h2
    pval = binomtest(h1, total_phased, p=0.5).pvalue if total_phased > 0 else None

    row = {
        "sample": sample,
        "locus": gene,
        "TOTAL_H1": h1,
        "TOTAL_H2": h2,
        "TOTAL_UNPHASED": unph,
        "UNANNOTATED_H1": unanno_h1,
        "UNANNOTATED_H2": unanno_h2,
        "UNANNOTATED_UNPHASED": unanno_unph,
        "BINOM_P": pval,
    }

    # Expand a few standard categories for convenience
    def add_pair(tag_prefix, pos_label, neg_label, alias_pos, alias_neg):
        row[f"H1_{alias_pos}"] = counts.get(("1", f"{tag_prefix}:{pos_label}"), 0)
        row[f"H1_{alias_neg}"] = counts.get(("1", f"{tag_prefix}:{neg_label}"), 0)
        row[f"H2_{alias_pos}"] = counts.get(("2", f"{tag_prefix}:{pos_label}"), 0)
        row[f"H2_{alias_neg}"] = counts.get(("2", f"{tag_prefix}:{neg_label}"), 0)
        row[f"UNPHASED_{alias_pos}"] = counts.get(("unphased", f"{tag_prefix}:{pos_label}"), 0)
        row[f"UNPHASED_{alias_neg}"] = counts.get(("unphased", f"{tag_prefix}:{neg_label}"), 0)

    # IN: in_annotation vs not_in_annotation
    add_pair("IN", "in_annotation", "not_in_annotation", "IN_ANNOTATION", "NOT_IN_ANNOTATION")
    # CD: coding vs non_coding
    add_pair("CD", "coding", "non_coding", "CODING", "NONCODING")
    # NM: NMD_sensitive vs not_NMD_sensitive
    add_pair("NM", "NMD_sensitive", "not_NMD_sensitive", "NMD", "NOT_NMD")
    # OR: has_ORF vs no_ORF
    add_pair("OR", "has_ORF", "no_ORF", "HAS_ORF", "NO_ORF")

    return row

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--tagged-bam', required=True)
    parser.add_argument('--passes', required=True)
    parser.add_argument('--bed', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--out', required=True)

    args = parser.parse_args()
    pass_genes = read_pass_genes(args.passes, args.sample, args.out)
    bed = read_bed_grouped(args.bed)
    bam = pysam.AlignmentFile(args.tagged_bam, "rb")


    rows = []
    # iterate only genes that PASS for this sample and exist in bed
    for gene in sorted(pass_genes & set(bed.keys())):
        gene_counts = defaultdict(int)
        for chrom, ivals in bed[gene].items():
            cc = count_in_gene(bam, chrom, ivals)
            for k, v in cc.items():
                gene_counts[k] += v
        rows.append(summarise_counts(args.sample, gene, gene_counts))

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
