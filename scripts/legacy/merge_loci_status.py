#!/usr/bin/env python3
"""
merge_loci_status.py
--------------------
Merge per-sample *.loci.csv files from IDENTIFY_POTENTIAL_LOCI into wide tables.

Input loci CSV columns expected (at minimum):
  sample, chrom, start, end, gene, ensembl_id, status
(Other columns may exist; ignored.)

Outputs:
  - successful_loci.all.csv  : all unique loci across all input files
  - successful_loci.pass.csv : only loci with >=1 PASS across samples

Definition of a locus key:
  (chrom, start, end, gene, ensembl_id)
We DO NOT use BED here. We trust the per-sample loci outputs.

Sample column names:
  derived from filename basename with .loci.csv/.csv stripped.

Missing locus status for a sample:
  MISSING (configurable below)
"""

import argparse
import csv
import os
from collections import defaultdict

MISSING = "MISSING"

KEY_FIELDS = ["chrom", "start", "end", "gene", "ensembl_id"]

def sample_name_from_path(p: str) -> str:
    base = os.path.basename(p)
    if base.endswith(".loci.csv"):
        return base[:-len(".loci.csv")]
    if base.endswith(".csv"):
        return base[:-len(".csv")]
    return base

def norm_row_keys(row: dict) -> dict:
    # normalise header keys like "status\r" or "Status " etc
    return { (k or "").strip().lower(): (v or "") for k, v in row.items() }

def locus_key(row: dict):
    # Keep as strings so we don’t accidentally change formatting
    return tuple((row.get(f, "") or "").strip() for f in KEY_FIELDS)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-all", default="successful_loci.all.csv", help="Wide table of all loci")
    ap.add_argument("--out-pass", default="successful_loci.pass.csv", help="Wide table filtered to loci with >=1 PASS")
    ap.add_argument("--missing", default=MISSING, help="Fill value when a locus is absent for a sample")
    ap.add_argument("csvs", nargs="+", help="Input per-sample loci CSVs (*.loci.csv)")
    args = ap.parse_args()

    # locus -> {sample -> status}
    status_map = defaultdict(dict)
    # locus -> canonical metadata (first-seen values)
    meta_map = {}

    samples = []
    for path in args.csvs:
        sname = sample_name_from_path(path)
        samples.append(sname)

        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                continue

            for raw in reader:
                row = norm_row_keys(raw)

                # Pull required fields
                # Note: your input has 'status' as the last col
                k = locus_key(row)
                if not any(k):
                    continue

                st = (row.get("status", "") or "").strip()
                st = st.upper() if st else ""

                if k not in meta_map:
                    meta_map[k] = dict(zip(KEY_FIELDS, k))

                status_map[k][sname] = st if st else args.missing

    # Unique samples (stable order)
    samples = sorted(set(samples))

    # Sort loci by genomic order where possible (chrom as string, start numeric-ish)
    def sort_key(k):
        chrom, start, end, gene, ens = k
        try:
            s = int(start)
        except:
            s = 0
        try:
            e = int(end)
        except:
            e = 0
        return (chrom, s, e, gene, ens)

    loci_sorted = sorted(meta_map.keys(), key=sort_key)

    header = KEY_FIELDS + samples

    # Build rows
    all_rows = []
    pass_rows = []

    for lk in loci_sorted:
        base = [meta_map[lk][f] for f in KEY_FIELDS]
        sts = [status_map[lk].get(s, args.missing) for s in samples]
        row = base + sts
        all_rows.append(row)

        if any(x == "PASS" for x in sts):
            pass_rows.append(row)

    # Write outputs
    with open(args.out_all, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow(header)
        w.writerows(all_rows)

    with open(args.out_pass, "w", newline="") as outfh:
        w = csv.writer(outfh)
        w.writerow(header)
        w.writerows(pass_rows)

    print(f"[INFO] Samples: {len(samples)}")
    print(f"[INFO] Loci total: {len(all_rows)}")
    print(f"[INFO] Loci with >=1 PASS: {len(pass_rows)}")
    print(f"[INFO] Wrote: {args.out_all}")
    print(f"[INFO] Wrote: {args.out_pass}")

if __name__ == "__main__":
    main()
