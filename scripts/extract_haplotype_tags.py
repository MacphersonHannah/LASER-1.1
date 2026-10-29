#!/usr/bin/env python3
import argparse, pysam, csv

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bams", nargs="+", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    BATCH = 100_000
    with open(args.out, "w", buffering=1024*1024) as out_f:
        out_f.write("read_id\thaplotype\n")
        batch = []
        for bam in args.bams:
            with pysam.AlignmentFile(bam, "rb") as bf:
                for read in bf.fetch(until_eof=True):
                    hp = None
                    if read.has_tag("HP"):
                        try:
                            hp = str(read.get_tag("HP"))
                        except Exception:
                            hp = None
                    batch.append(f"{read.query_name}\t{hp or 'unphased'}\n")
                    if len(batch) >= BATCH:
                        out_f.writelines(batch)
                        batch.clear()
        if batch:
            out_f.writelines(batch)
    print(f"[+] Wrote haplotype map to {args.out}")

if __name__ == "__main__":
    main()
