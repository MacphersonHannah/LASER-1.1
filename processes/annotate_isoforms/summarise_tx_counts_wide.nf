nextflow.enable.dsl = 2

// Not currently used. Legacy, kept for reference.
// output tsv columns: transcript_id, associated_gene_sqanti, total, H1, H2, UA

process COUNT_TX_PER_SAMPLE {
  tag "${sample}"

  publishDir "${params.outdir}/annotate_isoforms", mode: 'link'

  input:
    tuple val(sample), path(bam), path(bai)
    path lookup_sqlite

  output:
    tuple val(sample), path("${sample}.counts.tsv")

  errorStrategy 'terminate'
  maxRetries 1

  script:
  """
  set -euo pipefail

  # Extract read_id and HP tag from the tagged BAM
  # HP may be absent; mark those as UA ie unassigned
  samtools view -@ ${task.cpus} "${bam}" \
    | awk 'BEGIN{FS="\\t"}{
        hp="UA";
        for(i=12;i<=NF;i++){
          if(\$i ~ /^HP:i:/){
            split(\$i,a,":"); hp=a[3]
          }
        }
        print \$1"\\t"hp
      }' \
    > reads_hp.tsv

  # If read names in the BAM carry suffixes (eg /ccs) but SQLite uses bare IDs,
  # uncomment the next line to strip everything after the first space or slash so it matches properly
  # awk -F "\\t" '{gsub(/[\\s/].*/, "", \$1); print}' reads_hp.tsv > reads_hp.norm.tsv && mv reads_hp.norm.tsv reads_hp.tsv

  # Join to isoform_lookup in the SQLite to get transcript_id + associated_gene_sqanti,
  # then count per transcript_id split by haplotype.

  sqlite3 -cmd ".mode tabs" "${lookup_sqlite}" '
    PRAGMA temp_store=MEMORY;
    PRAGMA synchronous=OFF;
    PRAGMA journal_mode=OFF;
    DROP TABLE IF EXISTS reads;
    CREATE TEMP TABLE reads(read_id TEXT, hp TEXT);
    .import reads_hp.tsv reads

    WITH base AS (
      SELECT il.transcript_id,
             il.associated_gene_sqanti,
             COALESCE(NULLIF(r.hp, ""), "UA") AS hp
      FROM isoform_lookup il
      JOIN reads r USING(read_id)
    ),
    agg AS (
      SELECT transcript_id,
             associated_gene_sqanti,
             COUNT(*)                                        AS total,
             SUM(CASE WHEN hp = "1"  THEN 1 ELSE 0 END)      AS H1,
             SUM(CASE WHEN hp = "2"  THEN 1 ELSE 0 END)      AS H2,
             SUM(CASE WHEN hp = "UA" THEN 1 ELSE 0 END)      AS UA
      FROM base
      GROUP BY transcript_id, associated_gene_sqanti
    )
    SELECT transcript_id, associated_gene_sqanti, total, H1, H2, UA FROM agg
    ORDER BY transcript_id, associated_gene_sqanti;
  ' > "${sample}.counts.tsv"
  """
}

// Output is tx_counts_wide.tsv with columns: transcript_id, associated_gene_sqanti, <barcode>_total, <barcode>_H1, <barcode>_H2, <barcode>_UA  (for each barcode)

process MERGE_TX_COUNTS_WIDE {
  tag "merge-all-barcodes"
  publishDir "${params.outdir}/annotate_isoforms", mode: 'link', overwrite: true

  input:
    val tuples_list   // List of [sample, tsv path]
  
  output:
    path "tx_counts_wide.tsv"

  errorStrategy 'terminate'
  maxRetries 1

  script:
  """
  set -euo pipefail

  python - << 'PY'
import csv, sys, os
from collections import OrderedDict, defaultdict


PY
  # re-emit tuples_list into a tab file so Python can read it
  cat > tuples_list.tsv << 'EOF_TUP'
${tuples_list.collect{ it[0] + "\t" + it[1].toString() }.join("\n")}
EOF_TUP

  python - << 'PY'
import csv, sys, os
from collections import OrderedDict, defaultdict

# read the sample -> file mapping
pairs = []
with open("tuples_list.tsv") as f:
    for line in f:
        line = line.strip()
        if not line: 
            continue
        s, p = line.split("\\t", 1)
        pairs.append((s, p))

# key = (transcript_id, associated_gene_sqanti)
matrix = OrderedDict()
samples = [s for s,_ in pairs]
for s, path in pairs:
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\\t")
        for row in rdr:
            key = (row["transcript_id"], row["associated_gene_sqanti"])
            if key not in matrix:
                matrix[key] = {}
            # Coerce ints with fallback 0
            def to_int(x): 
                try: return int(x)
                except: return 0
            matrix[key][f"{s}_total"] = to_int(row.get("total", 0))
            matrix[key][f"{s}_H1"]    = to_int(row.get("H1", 0))
            matrix[key][f"{s}_H2"]    = to_int(row.get("H2", 0))
            matrix[key][f"{s}_UA"]    = to_int(row.get("UA", 0))

# build header
cols = ["transcript_id", "associated_gene_sqanti"]
for s in samples:
    cols.extend([f"{s}_total", f"{s}_H1", f"{s}_H2", f"{s}_UA"])

# write wide table
with open("tx_counts_wide.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\\t")
    w.writerow(cols)
    for (tx, gene), vals in matrix.items():
        row = [tx, gene]
        for s in samples:
            row.extend([
                vals.get(f"{s}_total", 0),
                vals.get(f"{s}_H1",    0),
                vals.get(f"{s}_H2",    0),
                vals.get(f"{s}_UA",    0),
            ])
        w.writerow(row)
PY
  """
}
