# inspired by a similar script from APTARS https://github.com/sid-sethi/APTARS
# last updated 22/07/2025

#!/usr/bin/env python

import sys
from Bio import SeqIO
import re
import json

primer_ids = []
five_p = ""
primers = sys.argv[1]

records = list(SeqIO.parse(primers, "fasta"))
for seq in records:
    if "_5p" in seq.id:
        five_p = seq.id.replace("-", "_")
    else:
        primer_ids.append(seq.id.replace("-", "_"))

if not five_p:
    raise ValueError("No primer ID with '_5p' found in primers FASTA. Required for demultiplexing naming.")

print(json.dumps({"five_p": five_p, "primer_ids": primer_ids}))
