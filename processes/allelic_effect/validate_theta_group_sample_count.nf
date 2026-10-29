process VALIDATE_THETA_GROUP_SAMPLE_COUNT {
  tag "theta_group_sample_count_check"
  conda "${projectDir}/envs_prebuilt/core_utils_env"
  publishDir "${params.outdir}/allelic_effect_testing/_validation", mode: 'link'
  input:
    val sample_list     // comma joined string of every sample name in this run
    path group_table    // CSV: sample,group

  output:
    path "theta_group_sample_count_ok.txt", emit: ok

  script:
  """
  set -euo pipefail
  python3 - "${sample_list}" "${group_table}" > theta_group_sample_count_ok.txt <<'PYEOF'
import sys, csv
from collections import defaultdict

samples = sys.argv[1].split(',')
group_table_path = sys.argv[2]

sample_to_group = {}
with open(group_table_path) as f:
    for row in csv.DictReader(f):
        sample_to_group[row['sample']] = row['group']

missing = [s for s in samples if s not in sample_to_group]
if missing:
    sys.exit(f"theta_group_table is missing group assignment(s) for: {', '.join(missing)}")

counts = defaultdict(int)
for s in samples:
    counts[sample_to_group[s]] += 1

for g, n in counts.items():
    if n < 2:
        sys.exit(f"Group '{g}' has only {n} sample(s) - calibrate_pooled_by_group needs at least 2 samples per group. Use 'lift' or 'fixed' for that group instead, or merge it into another group.")
    if n == 2:
        print(f"WARNING: group '{g}' has only 2 samples - calibrated theta rests on a very small basis.", file=sys.stderr)

print(f"{len(samples)} samples across {len(counts)} groups: " + ", ".join(f"{g}={n}" for g, n in counts.items()))
PYEOF
  """
}
