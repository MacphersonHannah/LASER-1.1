#!/usr/bin/env python3
# updated for speed 06/10/2025
# adapted from original script by 'etseng@pacb.com' in cDNA Cupcake

"""
After running collapse, combine with the cluster FL/nFL assignment to get abundance information.

cluster_report.csv format in IsoSeq3:

cluster_id,read_id,read_type
transcript/0,m54056_171130_193019/8388911/ccs,FL

where .group.txt is:

PB.1.1     transcript/0,transcript1
PB.1.2     transcript/2
"""
import os
import sys
import csv
import sqlite3
import tempfile
from collections import Counter
from argparse import ArgumentParser


def read_group_filename(group_filename):
    """
    Parse the *.group.txt (PB -> transcripts) into a dict:
      transcript_id (eg 'transcript/0') -> pbid (eg 'PB.1.1')
    """
    cid_info = {}
    with open(group_filename) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            pbid, members = line.split('\t', 1)
            # members may be "transcript/0,transcript/29"
            for cid in members.split(','):
                cid_info[cid] = pbid
    return cid_info


def sorted_pbids(pbid_iterable):
    """
    Sort PB IDs as before by numeric chunks after 'PB.' (PB.1.2 < PB.2.1).
    Example key: PB.10.2 -> [10, 2]
    """
    def key_fn(pbid):
        try:
            parts = pbid.split('.')[1:]  # drop 'PB'
            return [int(x) for x in parts]
        except Exception:
            return [pbid]
    return sorted(pbid_iterable, key=key_fn)


def write_abundance(abundance_path, tally_counter, write_header_comments=True):
    """
    Write abundance file in the same format:
      header comments +
      columns: pbid, count_fl, norm_fl
    Only counts FL unique reads (already tallied in 'tally_counter').
    """
    use_total_fl = sum(tally_counter.values())

    with open(abundance_path, 'w') as f:
        if write_header_comments:
            f.write("#\n")
            f.write("# -----------------\n")
            f.write("# Field explanation\n")
            f.write("# -----------------\n")
            f.write("# count_fl: Number of associated FL reads\n")
            f.write("# norm_fl: count_fl / total number of FL reads, mapped or unmapped\n")
            f.write(f"# Total Number of FL reads: {use_total_fl}\n")
            f.write("#\n")

        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(["pbid", "count_fl", "norm_fl"])

        for pbid in sorted_pbids(tally_counter.keys()):
            count_fl = tally_counter[pbid]
            norm_fl = (count_fl / use_total_fl) if use_total_fl else 0.0
            writer.writerow([pbid, count_fl, f"{norm_fl:.4e}"])


def write_mapped_fl_count(abundance_file, output_file):
    with open(abundance_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#') or line.startswith('pbid\t'):
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")




def process_cluster_report_fast(
    group_map,
    cluster_csv_path,
    read_stat_path,
    tally_counter,
    tmp_sqlite_path
):
    """
    One pass over cluster_report.csv, identical logic to the original:

    - FL + mapped transcript -> write unique/Y + tally pbid
    - FL + unmapped transcript -> write unmapped/Y
    - NonFL + mapped transcript -> store (read_id, pbid) to SQL (for later unique/ambiguous determination)
    - NonFL + unmapped transcript -> record read_id as potentially unmapped in SQL

    After the pass, we:
      - emit NonFL unique/ambiguous lines from SQL
      - emit NonFL unmapped (in SQL 'nonfl_unmapped' but not in 'nonfl_map') lines
    """
    # SQLite for NonFL aggregation (disk-backed, scales to huge inputs)
    con = sqlite3.connect(tmp_sqlite_path)
    cur = con.cursor()
    cur.execute("PRAGMA journal_mode=OFF")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA temp_store=MEMORY")
    cur.execute("PRAGMA cache_size=-500000")

    cur.execute("""
        CREATE TABLE IF NOT EXISTS nonfl_map (
            read_id TEXT,
            pbid    TEXT,
            PRIMARY KEY (read_id, pbid)
        )
    """)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS nonfl_unmapped (
            read_id TEXT PRIMARY KEY
        )
    """)
    con.commit()

    # Buffered writer for read_stat
    out_f = open(read_stat_path, 'w', buffering=1024*1024)
    out_f.write("id\tis_fl\tstat\tpbid\n")

    # Read header -> indexes
    with open(cluster_csv_path, newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        # Find required columns
        try:
            idx_cluster = header.index("cluster_id")
            idx_readid  = header.index("read_id")
            idx_readtyp = header.index("read_type")
        except ValueError as e:
            raise SystemExit(f"{cluster_csv_path}: missing required columns: {e}")

        # Batch buffers
        sql_batch_map = []
        sql_batch_unmapped = []
        WRITE_FLUSH = 100_000
        SQL_FLUSH = 100_000
        wbuf = []

        for row in reader:
            # robust guard for ragged lines
            if len(row) <= max(idx_cluster, idx_readid, idx_readtyp):
                continue

            cid = row[idx_cluster]
            rid = row[idx_readid]
            rtp = row[idx_readtyp]

            pbid = group_map.get(cid, None)

            if rtp == 'FL':
                if pbid is not None:
                    # mapped FL -> unique
                    wbuf.append(f"{rid}\tY\tunique\t{pbid}\n")
                    tally_counter[pbid] += 1
                else:
                    # unmapped FL
                    wbuf.append(f"{rid}\tY\tunmapped\tNA\n")
            else:
                # NonFL
                if pbid is not None:
                    sql_batch_map.append((rid, pbid))
                else:
                    sql_batch_unmapped.append((rid,))

            # periodic flush for speed/memory
            if len(wbuf) >= WRITE_FLUSH:
                out_f.writelines(wbuf)
                wbuf.clear()

            if len(sql_batch_map) >= SQL_FLUSH:
                cur.executemany(
                    "INSERT OR IGNORE INTO nonfl_map(read_id, pbid) VALUES (?,?)",
                    sql_batch_map
                )
                sql_batch_map.clear()

            if len(sql_batch_unmapped) >= SQL_FLUSH:
                cur.executemany(
                    "INSERT OR IGNORE INTO nonfl_unmapped(read_id) VALUES (?)",
                    sql_batch_unmapped
                )
                sql_batch_unmapped.clear()

        # final flush
        if wbuf:
            out_f.writelines(wbuf)
        if sql_batch_map:
            cur.executemany(
                "INSERT OR IGNORE INTO nonfl_map(read_id, pbid) VALUES (?,?)",
                sql_batch_map
            )
        if sql_batch_unmapped:
            cur.executemany(
                "INSERT OR IGNORE INTO nonfl_unmapped(read_id) VALUES (?)",
                sql_batch_unmapped
            )
        con.commit()

    # Emit NonFL unique/ambiguous from SQL (grouped by read_id)
    # Using GROUP_CONCAT is compact and fast
    for rid, concat_pbids in cur.execute(
        "SELECT read_id, GROUP_CONCAT(pbid) FROM nonfl_map GROUP BY read_id"
    ):
        pbids = concat_pbids.split(',') if concat_pbids else []
        if not pbids:
            continue
        stat = 'unique' if len(set(pbids)) == 1 else 'ambiguous'
        # Write one line per pbid (original behavior)
        for p in pbids:
            out_f.write(f"{rid}\tN\t{stat}\t{p}\n")

    # Emit NonFL unmapped (present in nonfl_unmapped but not in mapping)
    # LEFT JOIN to find those with NULL pbid
    for (rid,) in cur.execute("""
        SELECT u.read_id
        FROM nonfl_unmapped u
        LEFT JOIN nonfl_map m ON u.read_id = m.read_id
        WHERE m.read_id IS NULL
    """):
        out_f.write(f"{rid}\tN\tunmapped\tNA\n")

    out_f.close()
    con.close()


def get_abundance_post_collapse(collapse_prefix, cluster_report_csv, output_prefix):
    """
    # Wrapper matching the original function's behaviour and outputs.
    """
    group_filename = collapse_prefix + ".group.txt"
    if not os.path.exists(group_filename):
        print(f"File {group_filename} does not exist. Abort!", file=sys.stderr)
        sys.exit(-1)
    if not os.path.exists(cluster_report_csv):
        print(f"File {cluster_report_csv} does not exist. Abort!", file=sys.stderr)
        sys.exit(-1)

    cid_info = read_group_filename(group_filename)

    # Output paths (same as before)
    read_stat_path   = output_prefix + '.read_stat.txt'
    abundance_path   = output_prefix + '.abundance.txt'
    mapped_out_path  = output_prefix.replace('.collapsed', '_mapped_fl_count.txt')

    # Tally FL uniques while streaming cluster_report
    tally = Counter()

    # SQLite temp DB for NonFL aggregation (disk-backed)
    with tempfile.TemporaryDirectory() as tmpd:
        tmp_sqlite = os.path.join(tmpd, "nonfl.sqlite")
        process_cluster_report_fast(
            cid_info,
            cluster_report_csv,
            read_stat_path,
            tally,
            tmp_sqlite
        )

    # Abundance from tally (identical calculation as original)
    write_abundance(abundance_path, tally, write_header_comments=True)

    # Mapped FL count
    write_mapped_fl_count(abundance_path, mapped_out_path)

    # stderr messages (kept for parity)
    print("Read stat file written to", read_stat_path, file=sys.stderr)
    print("Abundance file written to", abundance_path, file=sys.stderr)
    print("Mapped FL count file written to", mapped_out_path, file=sys.stderr)

    return {
        "read_stat": read_stat_path,
        "abundance": abundance_path,
        "mapped": mapped_out_path
    }



if __name__ == "__main__":
    parser = ArgumentParser(
        "Get abundance/read stat information after running collapse script. Works for Iso-Seq1, 2, and 3 output."
    )
    parser.add_argument("collapse_prefix", help="Collapse prefix (must have .group.txt)")
    parser.add_argument("cluster_report", help="Cluster CSV report")
    args = parser.parse_args()

    sample = os.path.basename(args.collapse_prefix).replace('.collapsed', '').replace('_clustered', '')
    outdir = os.path.join("Collapsed_isoforms", sample)
    os.makedirs(outdir, exist_ok=True)
    output_prefix = os.path.join(outdir, sample + ".collapsed")

    result = get_abundance_post_collapse(args.collapse_prefix, args.cluster_report, output_prefix)

    # Print paths for Nextflow to parse (same keys as before)
    for name, path in result.items():
        print(f"{name}:{path}")