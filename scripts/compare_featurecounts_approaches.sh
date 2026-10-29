#!/usr/bin/env bash
# compare_featurecounts_approaches.sh
# for testing - Tests multiple featureCounts strategies for assigning PacBio FLNC
# reads to isocall-called transcripts.
set -euo pipefail

BAM=""
GTF=""
OUTDIR="fc_comparison"
THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)     BAM="$2";     shift 2 ;;
        --gtf)     GTF="$2";     shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$BAM" || -z "$GTF" ]] && { echo "ERROR: --bam and --gtf required"; exit 1; }
BAM_BASENAME=$(basename "$BAM")
mkdir -p "$OUTDIR"

echo "=== featureCounts transcript-assignment comparison ==="
echo "BAM:     $BAM"
echo "GTF:     $GTF"
echo "Outdir:  $OUTDIR"
echo ""


# Full summary

summarise_core() {
    local core_file="$1"
    local label="$2"
    awk -F'\t' -v label="$label" '
        /^#/ { next }
        {
            total++
            s = $2
            status[s]++

            # For ambiguous reads, bin by how many transcripts they are split between.
            # Column 3 (n_targets) is the count; column 4 lists them comma-separated.
            # featureCounts sets status="Ambiguous" when a read overlaps >1 feature
            # equally and cannot be resolved by largestOverlap.
            if (s == "Ambiguous") {
                n = $3 + 0
                if (n <= 1) n = 2   # guard: Ambiguous implies >=2
                if      (n == 2) amb[2]++
                else if (n == 3) amb[3]++
                else if (n == 4) amb[4]++
                else             amb[5]++   # 5 = "5 or more"
            }
        }
        END {
            assigned   = status["Assigned"]   + 0
            ambiguous  = status["Ambiguous"]  + 0
            nofeat     = status["Unassigned_NoFeatures"] + 0
            other      = total - assigned - ambiguous - nofeat

            printf "\n--- %s ---\n", label
            printf "  Total reads    : %d\n", total
            printf "  Assigned       : %d  (%.1f%%)\n", assigned,  100*assigned/total
            printf "  Ambiguous      : %d  (%.1f%%)\n", ambiguous, 100*ambiguous/total
            printf "  Dropped (NoFeat): %d  (%.1f%%)\n", nofeat,   100*nofeat/total
            printf "  Other unassigned: %d  (%.1f%%)\n", other,    100*other/total

            if (ambiguous > 0) {
                printf "\n  Ambiguity breakdown (how many transcripts each ambiguous read is split between):\n"
                printf "    2-way  : %d  (%.1f%% of ambiguous)\n", amb[2]+0, 100*(amb[2]+0)/ambiguous
                printf "    3-way  : %d  (%.1f%% of ambiguous)\n", amb[3]+0, 100*(amb[3]+0)/ambiguous
                printf "    4-way  : %d  (%.1f%% of ambiguous)\n", amb[4]+0, 100*(amb[4]+0)/ambiguous
                printf "    5-way+ : %d  (%.1f%% of ambiguous)\n", amb[5]+0, 100*(amb[5]+0)/ambiguous
            }
        }
    ' "$core_file"
}


# Run featureCounts, move CORE output, summarise

run_fc() {
    local label="$1"
    local extra_args="$2"
    local prefix="$OUTDIR/${label}"
    echo "--- Running: $label ---"
    featureCounts \
        -L \
        -T "$THREADS" \
        -a "$GTF" \
        -t exon -g transcript_id \
        --primary \
        -s 1 \
        -R CORE \
        $extra_args \
        "$BAM" \
        -o "${prefix}.fc.txt" \
        2>"${prefix}.fc.log"

    local core="${BAM_BASENAME}.featureCounts"
    if [[ ! -f "$core" ]]; then
        echo "ERROR: CORE file $core not found after $label" >&2
        return 1
    fi
    mv "$core" "${prefix}.core.tsv"
    summarise_core "${prefix}.core.tsv" "$label"
    echo ""
}

# Approaches


# A: largestOverlap only
run_fc "A_largestOverlap" "--largestOverlap"

# B: largestOverlap + fracOverlapFeature 0.5
run_fc "B_fracOverlapFeature_0.5" "--largestOverlap --fracOverlapFeature 0.5"

# C: largestOverlap + fracOverlap 0.5
run_fc "C_fracOverlap_0.5" "--largestOverlap --fracOverlap 0.5"

# D: multiOverlap (-O) - diagnostic only
run_fc "D_multiOverlap_diagnostic" "--largestOverlap -O"


# One-line comparison summary

for label in A_largestOverlap B_fracOverlapFeature_0.5 C_fracOverlap_0.5 D_multiOverlap_diagnostic; do
    core="$OUTDIR/${label}.core.tsv"
    [[ -f "$core" ]] || continue
    awk -F'\t' -v label="$label" '
        /^#/ { next }
        { total++; s[$2]++ }
        END {
            printf "%-40s  assigned=%.1f%%  ambiguous=%.1f%%  dropped=%.1f%%\n",
                label,
                100*(s["Assigned"]+0)/total,
                100*(s["Ambiguous"]+0)/total,
                100*(s["Unassigned_NoFeatures"]+0)/total
        }
    ' "$core"
done


# Write Assigned read>transcript TSVs for whichever approach you pick

for label in A_largestOverlap B_fracOverlapFeature_0.5 C_fracOverlap_0.5 D_multiOverlap_diagnostic; do
    core="$OUTDIR/${label}.core.tsv"
    out="$OUTDIR/${label}.read2transcript.tsv"
    [[ -f "$core" ]] || continue
    awk -F'\t' '
        /^#/ { next }
        $2 == "Assigned" && NF >= 4 {
            split($4, ids, ",")
            print $1 "\t" ids[1]
        }
    ' "$core" > "$out"
    echo "  $label: $(wc -l < "$out") assigned read>transcript rows > $out"
done


# How many reads does each strict filter drop vs baseline A?
echo ""
echo "Reads assigned in A but lost by stricter filters:"
echo "  Lost by B (fracOverlapFeature 0.5): $(comm -23 \
    <(awk '{print $1}' "$OUTDIR/A_largestOverlap.read2transcript.tsv" | sort) \
    <(awk '{print $1}' "$OUTDIR/B_fracOverlapFeature_0.5.read2transcript.tsv" | sort) \
    | wc -l)"
echo "  Lost by C (fracOverlap 0.5):        $(comm -23 \
    <(awk '{print $1}' "$OUTDIR/A_largestOverlap.read2transcript.tsv" | sort) \
    <(awk '{print $1}' "$OUTDIR/C_fracOverlap_0.5.read2transcript.tsv" | sort) \
    | wc -l)"

echo ""
echo "Done. Full CORE files and read2transcript TSVs in $OUTDIR/"