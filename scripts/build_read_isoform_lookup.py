"""
Creates a SQLite database linking reads > isoforms > SQANTI > ORFannotate 
Reads from read_sample_map, then haplotype/featureCounts/transcript assignments are attached where available.

Very important that every read is linked to a transcript, hence careful choice of transcript calling tools. Not all provide this info

Platforms: pacbio (cluster_report.csv + group.txt), ont (IsoQuant), flair (read.map.txt)
"""

import argparse
import csv
import sqlite3
from typing import Dict, List, Tuple
import duckdb

# helpers

def read_cluster_report(path: str) -> List[Dict[str, str]]:
    """Return list of dicts for cluster_report.csv (read_id + cluster_id)."""
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split(",")
        if "read_id" not in header or "cluster_id" not in header:
            raise ValueError(f"{path} missing required columns read_id / cluster_id")
        for line in fh:
            toks = line.strip().split(",")
            if len(toks) < len(header):
                continue
            row = dict(zip(header, toks))
            rows.append(row)
    return rows


def read_group_txt(path: str) -> Dict[str, str]:
    """Return mapping from transcript name > pbid."""
    m = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            toks = line.split("\t")
            if len(toks) >= 2:
                pbid, members = toks[0], toks[1:]
                for t in members:
                    m[t] = pbid
    return m


def read_orf_summary(path: str) -> Dict[str, Dict[str, str]]:
    """Return dict {pbid: dict of ORF fields}."""
    d = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            pbid = row.get("pbid") or row.get("transcript_id") or row.get("PBID")
            if pbid:
                d[pbid] = row
    return d


def read_sqanti_classification(path: str) -> Dict[str, Dict[str, str]]:
    """Return dict {pbid: dict of SQANTI fields}."""
    d = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            pbid = row.get("isoform") or row.get("pbid") or row.get("PBID")
            if pbid:
                d[pbid] = row
    return d


def read_sample_map(path: str) -> Dict[str, str]:
    """Return dict {read_id: sample_name}. This should define the read universe."""
    d = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rid = row.get("read_id") or row.get("ZMW") or row.get("name")
            smp = row.get("sample") or row.get("barcode") or row.get("Sample")
            if rid and smp:
                d[rid] = smp
    return d


def read_haplotypes(tsv_paths: List[str]) -> Dict[str, str]:
    """Merge multiple TSVs with read_id + haplotype columns."""
    h = {}
    for p in (tsv_paths or []):
        with open(p) as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if r.fieldnames and ('read_id' in r.fieldnames and 'haplotype' in r.fieldnames):
                for row in r:
                    rid = row.get('read_id')
                    hp = row.get('haplotype')
                    if rid:
                        h[rid] = hp
            else:
                fh.seek(0)
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    toks = line.strip().split("\t")
                    if len(toks) >= 2:
                        h[toks[0]] = toks[1]
    return h


def _status_rank(status: str) -> int:
    """
    Higher is better. We prefer Assigned over anything else.
    """
    if status == "Assigned":
        return 100
    if status.startswith("Unassigned_"):
        return 10
    if status in ("", "NA", None):
        return 0
    return 1


def read_featurecounts(tsv_paths: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Collapse featureCounts CORE TSV(s) to one row per read_id.
    Prefers Assigned rows; falls back to best Unassigned status.
    Returns dict: read_id -> dict(status, n_targets, genes)
    """
    tmp: Dict[str, List[Dict[str, str]]] = {}

    for p in (tsv_paths or []):
        with open(p) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                toks = line.rstrip("\n").split("\t")
                if len(toks) < 2:
                    continue

                rid = toks[0]
                status = toks[1] if len(toks) >= 2 else ""
                n_targets = toks[2] if len(toks) >= 3 else ""
                genes = toks[3] if len(toks) >= 4 else ""

                tmp.setdefault(rid, []).append({
                    "status": status,
                    "n_targets": n_targets,
                    "genes": genes
                })

    collapsed: Dict[str, Dict[str, str]] = {}
    for rid, recs in tmp.items():
        assigned = [r for r in recs if r.get("status") == "Assigned"]

        if assigned:
            gene_set = set()
            for r in assigned:
                gs = (r.get("genes") or "").strip()
                if gs and gs != "NA":
                    for g in gs.split(","):
                        g = g.strip()
                        if g:
                            gene_set.add(g)

            genes_out = ",".join(sorted(gene_set)) if gene_set else "NA"
            collapsed[rid] = {
                "status": "Assigned",
                "n_targets": str(len(gene_set)) if gene_set else "1",
                "genes": genes_out
            }
        else:
            best = max(recs, key=lambda r: _status_rank(r.get("status", "")))
            collapsed[rid] = {
                "status": best.get("status", ""),
                "n_targets": best.get("n_targets", ""),
                "genes": best.get("genes", "")
            }

    return collapsed


def build_pbid_by_read(cluster_rows: List[Dict[str, str]], trans_to_pbid: Dict[str, str]) -> Dict[str, Tuple[str, str]]:
    """
    Build mapping: read_id -> (transcript_name, pbid)
    cluster_report: read_id -> cluster_id (transcript_name)
    group_txt: transcript_name -> pbid
    """
    out: Dict[str, Tuple[str, str]] = {}
    for row in cluster_rows:
        rid = row.get("read_id")
        transcript = row.get("cluster_id", "")
        if not rid:
            continue
        pbid = trans_to_pbid.get(transcript, "")
        out[rid] = (transcript, pbid)
    return out


# SQLite writing

def write_transcripts_table(conn: sqlite3.Connection, transcripts: List[Dict[str, str]]) -> None:
    if not transcripts:
        # still create empty table with at least pbid
        cur = conn.cursor()
        cur.execute('CREATE TABLE transcripts ("pbid" TEXT)')
        cur.execute("CREATE UNIQUE INDEX idx_transcripts_pbid ON transcripts(pbid)")
        conn.commit()
        return

    cur = conn.cursor()
    cols = ", ".join(f'"{c}" TEXT' for c in transcripts[0].keys())
    cur.execute(f"CREATE TABLE transcripts ({cols})")
    cur.executemany(
        f"INSERT INTO transcripts VALUES ({','.join(['?'] * len(transcripts[0]))})",
        [list(d.values()) for d in transcripts]
    )
    cur.execute("CREATE UNIQUE INDEX idx_transcripts_pbid ON transcripts(pbid)")
    conn.commit()


def write_reads_table(conn: sqlite3.Connection, reads: List[Dict[str, str]]) -> None:
    if not reads:
        raise ValueError("No reads to write (reads_rows is empty). Check read_sample_map input.")

    cur = conn.cursor()
    cols = ", ".join(f'"{c}" TEXT' for c in reads[0].keys())
    cur.execute(f"CREATE TABLE reads ({cols})")
    # Insert in batches to avoid building one giant executemany call
    BATCH = 100_000
    row_template = f"INSERT INTO reads VALUES ({','.join(['?'] * len(reads[0]))})"
    for i in range(0, len(reads), BATCH):
        cur.executemany(row_template, [list(d.values()) for d in reads[i:i+BATCH]])
        conn.commit()

    # indexes + uniqueness guarantees
    cur.execute("CREATE INDEX idx_reads_read_id ON reads(read_id)")
    cur.execute("CREATE INDEX idx_reads_sample ON reads(sample)")
    cur.execute("CREATE INDEX idx_reads_pbid ON reads(pbid)")
    cur.execute("CREATE UNIQUE INDEX idx_reads_sample_read ON reads(sample, read_id)")
    conn.commit()


def read_isoquant_model_reads(path: str):
    """
    Parse IsoQuant transcript_model_reads.tsv.gz for ONT.
    Returns cluster_rows [{read_id, cluster_id}] and trans_to_pbid {isoform_id: isoform_id}.
    Columns: read_id <tab> isoform_id <tab>
    """
    import gzip
    cluster_rows = []
    trans_to_pbid = {}
    with open(path, 'rb') as _f:
        magic = _f.read(2)
    is_gz = magic == b'\x1f\x8b'
    opener = gzip.open if is_gz else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            toks = line.strip().split('\t')
            if len(toks) < 2:
                continue
            read_id = toks[0]
            if read_id == 'read_id':  # skip header line
                continue
            isoform_id = toks[1]  # column 2 = transcript_id
            if isoform_id in ('.', '*'):   # skip unassigned reads
                continue
            cluster_rows.append({'read_id': read_id, 'cluster_id': isoform_id})
            trans_to_pbid[isoform_id] = isoform_id
    return cluster_rows, trans_to_pbid


def read_flair_map(path: str):
    """
    Parse FLAIR quantify isoform.read.map.txt (isoform_id <tab> read1,read2,...).
    Strips position prefix from read IDs to match BAM QNAMEs.
    Returns cluster_rows and trans_to_pbid.

    NB FLAIR read IDs have the format: alignstart:alignend|uuid
    eg 214:796|8ba0174d-e6c0-414c-9ea4-0a0bb498e687
    This is different to isoseq hence separate handling
    """
    cluster_rows = []
    trans_to_pbid = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            toks = line.strip().split('\t')
            if len(toks) < 2:
                continue
            isoform_id = toks[0]
            reads_csv  = toks[1]
            if isoform_id in ('.', '*', '') or reads_csv in ('.', '*', ''):
                continue
            # Bug fix june 2026: FLAIR quantify appends gene/locus suffix to all isoform IDs in the read map
            # butSQANTI receives isoform IDs from the GTF without the suffix so they can mismatch.
            # This strips the last underscore-delimited segment to match SQANTI.
            # Rule: strip last '_segment' from any isoform_id containing '_'
            if '_' in isoform_id:
                pbid = isoform_id.rsplit('_', 1)[0]
            else:
                pbid = isoform_id
            trans_to_pbid[isoform_id] = pbid
            for read_id in reads_csv.split(','):
                read_id = read_id.strip()
                if not read_id:
                    continue
                # Strip FLAIR position prefix: "alignstart:alignend|uuid" -> "uuid"
                if '|' in read_id:
                    read_id = read_id.split('|', 1)[1]
                cluster_rows.append({'read_id': read_id, 'cluster_id': isoform_id})
    return cluster_rows, trans_to_pbid


#  main logic

def main():
    ap = argparse.ArgumentParser(description="Build read>isoform>annotation SQLite lookup (read universe from read_sample_map)")
    ap.add_argument("--cluster_report", required=True)
    ap.add_argument("--group_txt", required=True)
    ap.add_argument("--orf_summary", required=True)
    ap.add_argument("--sqanti_classification", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--read_sample_map", required=True)
    ap.add_argument("--haplotypes", action="append", default=[], help="TSV(s) with read_id + haplotype")
    ap.add_argument("--out", required=True)
    ap.add_argument("--featurecounts", action="append", default=[], help="TSV(s) with read_id + status + n_targets + genes (featureCounts)")
    ap.add_argument("--platform", default="pacbio", choices=["pacbio", "ont", "flair"])
    ap.add_argument("--memory_limit", default="500GB")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--temp_dir", default=".")

    args = ap.parse_args()

    print(f"[INFO] Building lookup for sample: {args.sample}")
    print(f"[INFO] Platform: {args.platform}")
    print(f"[INFO] Haplotype files: {len(args.haplotypes)}")
    print(f"[INFO] featureCounts CORE files: {len(args.featurecounts)}")

    orf_by_pbid = read_orf_summary(args.orf_summary)
    sqanti_by_pbid = read_sqanti_classification(args.sqanti_classification)

    # transcripts table (PBID-only) 
    transcript_rows = {}
    for pbid, orf_row in orf_by_pbid.items():
        merged = {"pbid": pbid}
        for k, v in orf_row.items():
            merged[f"orf_{k}"] = v
        sq_row = sqanti_by_pbid.get(pbid, {})
        for k, v in sq_row.items():
            merged[f"sqanti_{k}"] = v
        transcript_rows[pbid] = merged

    # also include SQANTI-only PBIDs that have no ORF record
    for pbid, sq_row in sqanti_by_pbid.items():
        if pbid in transcript_rows:
            continue
        merged = {"pbid": pbid}
        for k, v in sq_row.items():
            merged[f"sqanti_{k}"] = v
        transcript_rows[pbid] = merged

    # Write SQLite (stream reads directly - never build full reads_rows list)
    print(f"[INFO] Writing SQLite database > {args.out}")
    conn = sqlite3.connect(args.out)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    # Reduced cache - we're streaming so don't need 200MB of cache
    conn.execute("PRAGMA cache_size=-32000")

    write_transcripts_table(conn, list(transcript_rows.values()))
    del transcript_rows  # free memory before writing reads

    # Reads are counted per (sample, pbid, haplotype) rather than stored one row
    # per read - extract_transcript_summary_from_sqlite.py only ever asked the
    # reads table for that GROUP BY. DuckDB does the read_id joins on disk so
    # peak memory is bounded by --memory_limit, not by the number of reads.
    rid = "regexp_replace(read_id, '^[^|]*\\|', '')" if args.platform == 'flair' else "read_id"

    con = duckdb.connect()
    con.execute(f"SET memory_limit='{args.memory_limit}'")
    con.execute(f"SET threads={args.threads}")
    con.execute(f"SET temp_directory='{args.temp_dir}'")
    con.execute("SET preserve_insertion_order=false")

    # read universe + sample, one row per read
    con.execute(f"""
        CREATE TEMP TABLE rsm AS
        SELECT {rid} AS read_id, min(sample) AS sample
        FROM read_csv(?, delim='\t', header=true, all_varchar=true,
                      columns={{'read_id':'VARCHAR','sample':'VARCHAR'}})
        WHERE read_id IS NOT NULL AND sample IS NOT NULL
        GROUP BY 1
    """, [args.read_sample_map])

    # read>pbid, same per-platform parsing as read_flair_map / read_isoquant_model_reads
    if args.platform == 'flair':
        con.execute(f"""
            CREATE TEMP TABLE iso AS
            SELECT regexp_replace(r, '^[^|]*\\|', '') AS read_id,
                   min(CASE WHEN contains(isoform_id, '_')
                            THEN regexp_replace(isoform_id, '_[^_]*$', '')
                            ELSE isoform_id END) AS pbid
            FROM (
              SELECT isoform_id, unnest(string_split(reads_csv, ',')) AS r
              FROM read_csv(?, delim='\t', header=false, all_varchar=true,
                            max_line_size=536870912,
                            columns={{'isoform_id':'VARCHAR','reads_csv':'VARCHAR'}})
              WHERE isoform_id NOT IN ('.', '*', '') AND isoform_id IS NOT NULL
                AND reads_csv  NOT IN ('.', '*', '') AND reads_csv  IS NOT NULL
            )
            WHERE r IS NOT NULL AND r <> ''
            GROUP BY 1
        """, [args.cluster_report])
    elif args.platform == 'ont':
        con.execute("""
            CREATE TEMP TABLE iso AS
            SELECT read_id, min(isoform_id) AS pbid
            FROM read_csv(?, delim='\t', header=true, all_varchar=true,
                          columns={'read_id':'VARCHAR','isoform_id':'VARCHAR'})
            WHERE read_id IS NOT NULL AND isoform_id NOT IN ('.', '*', '')
            GROUP BY 1
        """, [args.cluster_report])
    else:
        # pacbio: same parsing as read_cluster_report() + read_group_txt().
        # cluster_report.csv is read by column name; group.txt is the pbid in field 1
        # and transcript names in the remaining tab-separated fields.
        con.execute("""
            CREATE TEMP TABLE trans_to_pbid AS
            SELECT tname, min(pbid) AS pbid
            FROM (
              SELECT trim(parts[1]) AS pbid, unnest(parts[2:]) AS tname
              FROM (
                SELECT string_split(trim(line, ' \t\r\n'), chr(9)) AS parts
                FROM read_csv(?, delim=chr(31), header=false, all_varchar=true,
                              quote='', escape='', columns={'line':'VARCHAR'})
                WHERE trim(line) <> ''
              )
              WHERE len(parts) >= 2
            )
            WHERE tname IS NOT NULL AND tname <> ''
            GROUP BY 1
        """, [args.group_txt])

        con.execute("""
            CREATE TEMP TABLE iso AS
            SELECT c.read_id, min(COALESCE(g.pbid, '')) AS pbid
            FROM read_csv(?, delim=',', header=true, all_varchar=true,
                          ignore_errors=true) c
            LEFT JOIN trans_to_pbid g ON g.tname = c.cluster_id
            WHERE c.read_id IS NOT NULL
            GROUP BY 1
        """, [args.cluster_report])

    if args.haplotypes:
        con.execute(f"""
            CREATE TEMP TABLE hap AS
            SELECT {rid} AS read_id, min(haplotype) AS hp
            FROM read_csv(?, delim='\t', header=true, all_varchar=true,
                          columns={{'read_id':'VARCHAR','haplotype':'VARCHAR'}})
            WHERE read_id IS NOT NULL
            GROUP BY 1
        """, [args.haplotypes])
    else:
        con.execute("CREATE TEMP TABLE hap AS SELECT NULL AS read_id, NULL AS hp WHERE false")

    rel = con.execute("""
        SELECT s.sample,
               COALESCE(i.pbid, '')         AS pbid,
               COALESCE(h.hp, 'not_in_bam') AS haplotype,
               COUNT(*)                     AS n
        FROM rsm s
        LEFT JOIN iso i USING (read_id)
        LEFT JOIN hap h USING (read_id)
        GROUP BY 1, 2, 3
    """)

    cur = conn.cursor()
    cur.execute('CREATE TABLE read_counts ("sample" TEXT, "pbid" TEXT, "haplotype" TEXT, "n" INTEGER)')
    n_rows = 0
    n_total = 0
    while True:
        batch = rel.fetchmany(100_000)
        if not batch:
            break
        cur.executemany("INSERT INTO read_counts VALUES (?,?,?,?)", batch)
        conn.commit()
        n_rows += len(batch)
        n_total += sum(r[3] for r in batch)

    cur.execute("CREATE UNIQUE INDEX idx_read_counts_key ON read_counts(sample, pbid, haplotype)")
    conn.commit()
    conn.close()
    con.close()

    print(f"[INFO] Counted {n_total:,} reads from read_sample_map (ALL reads)")
    print(f"[INFO] Wrote {n_rows:,} (sample, pbid, haplotype) rows to SQLite.")
    print("[INFO] Done.")

if __name__ == "__main__":
    main()
