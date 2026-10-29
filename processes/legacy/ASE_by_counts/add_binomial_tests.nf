process ADD_BINOMIAL_TESTS {
   publishDir "${params.outdir}/ase_by_counts/all", mode: 'copy', overwrite: true
  input:
    path merged_tsv

  output:
    path "ALL_COUNTS.with_tests.tsv", emit: with_tests

  // needs pandas + scipy in your ase env; same one you used for counts
  conda "${projectDir}/envs_prebuilt/ase_laser"

  script:
  """
  python3 - << 'PY'
import math
import pandas as pd

# SciPy compatibility across versions
try:
    from scipy.stats import binomtest as _binom
    def binom_p(h1, n): return _binom(h1, n, 0.5, alternative='two-sided').pvalue
except Exception:
    from scipy.stats import binom_test as _binom
    def binom_p(h1, n): return _binom(h1, n, 0.5)

df = pd.read_csv("${merged_tsv}", sep="\\t")

pairs = {
    "TOTAL": ("TOTAL_H1", "TOTAL_H2"),
    "IN_ANNOTATION": ("H1_IN_ANNOTATION", "H2_IN_ANNOTATION"),
    "NOT_IN_ANNOTATION": ("H1_NOT_IN_ANNOTATION", "H2_NOT_IN_ANNOTATION"),
    "CODING": ("H1_CODING", "H2_CODING"),
    "NONCODING": ("H1_NONCODING", "H2_NONCODING"),
    "NMD": ("H1_NMD", "H2_NMD"),
    "NOT_NMD": ("H1_NOT_NMD", "H2_NOT_NMD"),
    "HAS_ORF": ("H1_HAS_ORF", "H2_HAS_ORF"),
    "NO_ORF": ("H1_NO_ORF", "H2_NO_ORF"),
}

def compute_stats(row, a, b):
    h1 = int(row.get(a, 0) or 0)
    h2 = int(row.get(b, 0) or 0)
    n = h1 + h2
    if n == 0:
        return (None, None, None)
    try:
        p = float(binom_p(h1, n))
    except Exception:
        p = None
    af = h1 / n
    log2or = math.log2((h1 + 0.5) / (h2 + 0.5))  # pseudocount
    return (p, af, log2or)

# Run binomial tests for all categories
for cat, (a, b) in pairs.items():
    ps, afs, ors = [], [], []
    for _, r in df.iterrows():
        p, af, lr = compute_stats(r, a, b)
        ps.append(p); afs.append(af); ors.append(lr)
    df[f"BINOM_P_{cat}"] = ps
    df[f"AF_{cat}"]      = afs
    df[f"LOG2OR_{cat}"]  = ors

# Add ratings
def aste_rating(p):
    if p is None: return "no"
    if p < 0.01: return "significant"
    elif p < 0.03: return "almost"
    elif p < 0.1: return "trending"
    else: return "no"

df.insert(1, "asTE_rating", df["BINOM_P_TOTAL"].apply(aste_rating))
df.insert(2, "asTU_rating", "not_assessed")

# nice column order: original first, then tests grouped by category
test_cols = []
for cat in pairs.keys():
    test_cols += [f"BINOM_P_{cat}", f"AF_{cat}", f"LOG2OR_{cat}"]

# keep original order + tests appended
out_cols = list(df.columns.drop(test_cols, errors='ignore')) + test_cols
df[out_cols].to_csv("ALL_COUNTS.with_tests.tsv", sep="\\t", index=False)
PY
  """
}
