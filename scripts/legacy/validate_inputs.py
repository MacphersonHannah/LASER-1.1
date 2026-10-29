#!/usr/bin/env python3
import sys, os
from Bio import SeqIO

analysis_name, ccs_bam, primers = sys.argv[1:4]

ok = True
errors = []

# Check analysis_name
if " " in analysis_name:
    errors.append(f"ERROR: analysis_name '{analysis_name}' contains spaces.")

# Check BAM
if not os.path.exists(ccs_bam):
    errors.append(f"ERROR: BAM file '{ccs_bam}' does not exist.")
elif os.path.getsize(ccs_bam) == 0:
    errors.append(f"ERROR: BAM file '{ccs_bam}' is empty.")
elif not ccs_bam.endswith(".bam"):
    errors.append(f"ERROR: BAM file '{ccs_bam}' does not end with .bam extension.")

# Check primers
if not os.path.exists(primers):
    errors.append(f"ERROR: Primers file '{primers}' does not exist.")
else:
    bad_ids = []
    n_records = 0
    for rec in SeqIO.parse(primers, "fasta"):
        n_records += 1
        if not (rec.id.endswith("_5p") or rec.id.endswith("_3p")):
            bad_ids.append(rec.id)
    if n_records == 0:
        errors.append(f"ERROR: Primers file '{primers}' contains no records.")
    elif bad_ids:
        errors.append(f"ERROR: Primer IDs without _5p/_3p: {', '.join(bad_ids)}")

# Print errors or PASS
if errors:
    for err in errors:
        print(err)
    sys.exit(1)
else:
    sys.exit(0)