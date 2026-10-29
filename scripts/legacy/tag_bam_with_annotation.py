#!/usr/bin/env python3

### 06/10/2025 changes to speed up exponentially with large datasets

import argparse
import sqlite3
import pysam
from itertools import islice

CHUNK_SIZE = 20000  # adjust for memory vs speed; 20k-50k is good

def fetch_chunk(cur, read_ids):
    """Fetch rows for a batch of read IDs and return dict {read_id: row_dict}."""
    if not read_ids:
        return {}
    placeholders = ",".join("?" for _ in read_ids)
    sql = f"SELECT * FROM isoform_lookup WHERE read_id IN ({placeholders})"
    cur.execute(sql, tuple(read_ids))
    rows = cur.fetchall()
    cols = [d[0] for d in cur.description]
    return {r[cols.index("read_id")]: dict(zip(cols, r)) for r in rows}

def chunked_iterable(iterable, size):
    """Yield successive chunks of given size."""
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sqlite", required=True)
    ap.add_argument("--out_bam", required=True)
    args = ap.parse_args()

    conn = sqlite3.connect(args.sqlite)
    cur = conn.cursor()

    in_bam = pysam.AlignmentFile(args.bam, "rb")
    out_bam = pysam.AlignmentFile(args.out_bam, "wb", template=in_bam)

    # iterate in chunks
    reads_iter = in_bam.fetch(until_eof=True)
    for read_chunk in chunked_iterable(reads_iter, CHUNK_SIZE):
        # Get all read names for this batch
        read_ids = [r.query_name for r in read_chunk]
        lookup = fetch_chunk(cur, read_ids)

        for read in read_chunk:
            row = lookup.get(read.query_name)
            if row:
                in_ann = "in_annotation" if (row.get("reference_match_sqanti") or row.get("associated_transcript_sqanti")) else "not_in_annotation"
                read.set_tag("IN", in_ann, value_type="Z")

                cd = "coding" if row.get("coding_class") == "coding" else "non_coding"
                read.set_tag("CD", cd, value_type="Z")

                nm = "NMD_sensitive" if str(row.get("NMD_sensitive")).upper() == "TRUE" else "not_NMD_sensitive"
                read.set_tag("NM", nm, value_type="Z")

                has_orf = "has_ORF" if str(row.get("has_ORF")).upper() == "TRUE" else "no_ORF"
                read.set_tag("OR", has_orf, value_type="Z")
            else:
                for tag in ("IN", "CD", "NM", "OR"):
                    read.set_tag(tag, "read_unannotated", value_type="Z")

            out_bam.write(read)

    in_bam.close()
    out_bam.close()
    conn.close()
    pysam.index(args.out_bam)

if __name__ == "__main__":
    main()
