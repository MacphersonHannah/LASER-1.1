#!/usr/bin/env python3
"""
Identify isoforms that are haplotype-exclusive within a (barcode,locus)
according to three criteria:

    1. Present in exactly one haplotype (HP1 or HP2).  
    2. Read-count  >= --count-thresh   (default 3).  
    3. Within-haplotype proportion >= --prop-thresh (default 0.01).

Outputs one tab-separated file per barcode with the columns barcode  locus  isoform_id  haplotype  count  proportion
"""

import argparse
import csv
import sqlite3
from collections import defaultdict

 
# command-line parsing
 
def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db",          required=True, help="SQLite file built earlier")
    p.add_argument("--barcode",     required=True, help="Exact barcode string, eg BC01")
    p.add_argument("--loci",        required=True, nargs="+",
                   help="List of loci (gene names) that passed 'potential loci'")
    p.add_argument("--count-thresh", type=int, default=3)
    p.add_argument("--prop-thresh",  type=float, default=0.01)
    p.add_argument("--out",         required=True, help="Output TSV")
    return p.parse_args()


 
# helper queries
 
READS_Q = """
SELECT isoform_id,
       haplotype               AS hp,
       COUNT(*)                AS n
FROM   reads                   -- <- table created by BUILD_READ_ISOFORM_LOOKUP
WHERE  barcode   = ?
  AND  locus     = ?
  AND  hp IN ('1','2')         -- ignore unphased
GROUP BY isoform_id, hp
"""

TOTALS_Q = """
SELECT haplotype, COUNT(*) AS n_total
FROM   reads
WHERE  barcode = ?
  AND  locus   = ?
  AND  haplotype IN ('1','2')
GROUP BY haplotype
"""


 
# main work
 
def haplotype_unique(conn, barcode, locus, n_thresh, prop_thresh):
    """
    Yield (isoform_id, hp, count, proportion) rows that satisfy
    the 3-point definition for a single locus in one barcode.
    """
    cur = conn.cursor()

    # total reads per haplotype for this locus
    totals = {hp: n for hp, n in cur.execute(TOTALS_Q, (barcode, locus))}
    if not totals:
        return  # no phased reads

    # counts per isoform+haplotype
    rows = cur.execute(READS_Q, (barcode, locus))

    # bucket by isoform
    by_iso_hp = defaultdict(dict)        # isoform_id -> {hp: n}
    for iso, hp, n in rows:
        by_iso_hp[iso][hp] = n

    # apply the criteria
    for iso, counts in by_iso_hp.items():
        if len(counts) != 1:
            continue                     # present in both -> not exclusive

        hp, n = next(iter(counts.items()))

        if n < n_thresh:
            continue

        prop = n / totals[hp]
        if prop < prop_thresh:
            continue

        yield iso, hp, n, prop


def main():
    args = get_args()

    conn = sqlite3.connect(args.db)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["barcode", "locus", "isoform_id", "haplotype", "count", "proportion"]
        )

        for locus in args.loci:
            for iso, hp, n, prop in haplotype_unique(
                conn,
                args.barcode,
                locus,
                args.count_thresh,
                args.prop_thresh,
            ):
                w.writerow([args.barcode, locus, iso, hp, n, f"{prop:.5f}"])

    conn.close()


if __name__ == "__main__":
    main()
