// Per-sample quantification against the combined FLAIR transcriptome.
// Uses each sample's original reads so assignments are sample-specific.
// --stringent removed to avoid asymmetric allelic read loss in ASE.
// Direct RNA: see flair_quantify_drna.nf
//
// FLAIR_QUANTIFY is now a subworkflow with the same inputs and outputs as the old
// FLAIR_QUANTIFY process, so ont_cdna_preprocess.nf and
// ont_cdna_minimap_entry_preprocess.nf call it unchanged.
//
// params.flair_quant_chunks (default 24)
//   > 0 : parallel. `flair quantify` is a wrapper that pipes two commands
//         (flair_quantify.py):
//           minimap2 --MD -a -N 4 ... <isoforms.fa> <reads> | filter_transcriptome_align.py -s - ...
//         and the filter script reads every alignment in ONE process (~1 core per
//         sample). Here the same two commands run as separate Nextflow processes:
//           FLAIR_QUANT_INDEX  : minimap2 index of the isoforms, built once
//           FLAIR_QUANT_ALIGN  : the same minimap2 command, output dealt (whole reads)
//                                into N chunk BAMs
//           FLAIR_QUANT_FILTER : the same filter command, one task per chunk
//           FLAIR_QUANT_MERGE  : sum counts + join maps -> same two files, same folder
//         FLAIR assigns every read using only that read's own alignments, so the
//         result is identical to `flair quantify` (tested: byte-identical counts,
//         identical read->isoform assignments).
//     0 : FLAIR_QUANTIFY_SINGLE, the original `flair quantify` process, unchanged.
//
// INDEX: minimap2 splits its index when the reference exceeds -I (default 8G);
// a split index writes no @SQ lines and the filter fails. FLAIR_QUANT_INDEX uses
// a large -I so it is always one part.

nextflow.enable.dsl = 2


workflow FLAIR_QUANTIFY {

  take:
    reads          // (sample, fastq)
    isoforms_fa    // path, from FLAIR_TRANSCRIPTOME
    isoforms_bed   // path, from FLAIR_TRANSCRIPTOME

  main:
    // explicit null check: `?:` would turn a user-set 0 into the default
    def nchunks = (params.flair_quant_chunks != null) ? (params.flair_quant_chunks as int) : 24

    if (nchunks > 0) {
        log.info "FLAIR quantify: parallel, ${nchunks} chunks per sample"

        FLAIR_QUANT_INDEX(isoforms_fa)

        FLAIR_QUANT_ALIGN(reads, FLAIR_QUANT_INDEX.out.mmi, nchunks)

        // (sample, [chunk1..chunkN]) -> (sample, chunk) per chunk
        FLAIR_QUANT_ALIGN.out.chunks
            .flatMap { sample, chunks ->
                def list = (chunks instanceof List) ? chunks : [chunks]
                list.collect { c -> tuple(sample, c) }
            }
            .set { chunk_ch }

        FLAIR_QUANT_FILTER(chunk_ch, isoforms_bed)

        // merge a sample as soon as all its chunks are done; sorted so the merge
        // task hash is stable across -resume
        FLAIR_QUANT_FILTER.out.result
            .groupTuple(by: 0, size: nchunks, sort: true)
            .set { grouped_ch }

        FLAIR_QUANT_MERGE(grouped_ch)

        counts_ch   = FLAIR_QUANT_MERGE.out.counts
        read_map_ch = FLAIR_QUANT_MERGE.out.read_map

    } else {
        log.info "FLAIR quantify: original single process (flair quantify)"

        FLAIR_QUANTIFY_SINGLE(reads, isoforms_fa, isoforms_bed)

        counts_ch   = FLAIR_QUANTIFY_SINGLE.out.counts
        read_map_ch = FLAIR_QUANTIFY_SINGLE.out.read_map
    }

  emit:
    counts   = counts_ch
    read_map = read_map_ch
}


// ---------------------------------------------------------------------------
// ORIGINAL process (used when flair_quant_chunks = 0). Script unchanged.
// ---------------------------------------------------------------------------
process FLAIR_QUANTIFY_SINGLE {

    tag { sample }
    publishDir "${params.outdir}/preprocessing/flair/quantify/${sample}", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(fastq)   // per-sample reads (post-pychopper)
        path(isoforms_fa)                // from FLAIR_TRANSCRIPTOME
        path(isoforms_bed)               // from FLAIR_TRANSCRIPTOME

    output:
        tuple val(sample), path("${sample}.flair.counts.tsv"),     emit: counts
        tuple val(sample), path("${sample}.isoform.read.map.txt"), emit: read_map

    script:
    """
    # Manifest with fixed short condition/batch labels (no underscores).
    # FLAIR names the map file: output.sample_condition_batch.isoform.read.map.txt
    echo -e "${sample}\tcond\tb1\t${fastq}" > manifest.tsv

    flair quantify \\
        --isoforms ${isoforms_fa} \\
        --reads_manifest manifest.tsv \\
        --isoform_bed ${isoforms_bed} \\
        --output ${sample}.flair \\
        --sample_id_only \\
        --quality 0 \\
        --trust_ends \\
        --check_splice \\
        --generate_map \\
        --threads ${task.cpus}

    # Rename map file to predictable name - use wildcard in case FLAIR version
    # changes the exact suffix format
    mv ${sample}.flair.*.isoform.read.map.txt \\
       ${sample}.isoform.read.map.txt
    """
}


// ---------------------------------------------------------------------------
// PARALLEL processes
// ---------------------------------------------------------------------------

process FLAIR_QUANT_INDEX {

    tag { 'transcriptome' }

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        path(isoforms_fa)

    output:
        path("flair.isoforms.mmi"), emit: mmi

    script:
    def index_size = params.flair_quant_index_size ?: '100G'
    """
    set -euo pipefail
    # Same default k/w as flair quantify (no -x preset). -I only controls when the
    # index is split into parts; set high so it is always one part.
    minimap2 -I ${index_size} -t ${task.cpus} -d flair.isoforms.mmi ${isoforms_fa}
    """

    stub:
    """
    touch flair.isoforms.mmi
    """
}


process FLAIR_QUANT_ALIGN {

    tag { sample }

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(fastq)
        path(mmi)
        val(nchunks)

    output:
        tuple val(sample), path("chunks/${sample}.chunk_*.bam"), emit: chunks

    script:
    """
    set -euo pipefail
    mkdir -p chunks

    # Align once at full speed, then deal whole reads round-robin into ${nchunks} BAMs.
    # minimap2 writes all records of a read consecutively, and
    # no sorting is done, so each read's alignments stay together in one chunk.
    # Header lines go to every chunk.
    minimap2 --MD -a -N 4 -t ${task.cpus} ${mmi} ${fastq} \\
    | awk -v n=${nchunks} -v pre="chunks/${sample}.chunk_" '
        BEGIN { for (i = 0; i < n; i++) cmd[i] = sprintf("samtools view -b -o %s%03d.bam -", pre, i) }
        /^@/  { for (i = 0; i < n; i++) print | cmd[i]; next }
              { if (\$1 != last) { k = (k + 1) % n; last = \$1 } print | cmd[k] }
        END   { for (i = 0; i < n; i++) close(cmd[i]) }'

    # awk does not propagate samtools failures, so check every chunk explicitly
    n_found=\$(ls chunks/${sample}.chunk_*.bam | wc -l)
    if [ "\$n_found" -ne ${nchunks} ]; then
        echo "Expected ${nchunks} chunks, found \$n_found" >&2; exit 1
    fi
    samtools quickcheck -v chunks/${sample}.chunk_*.bam
    """

    stub:
    """
    mkdir -p chunks
    for i in \$(seq 0 \$((${nchunks} - 1))); do touch chunks/${sample}.chunk_\$(printf %03d \$i).bam; done
    """
}


process FLAIR_QUANT_FILTER {

    tag { chunk.baseName }

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(chunk)
        path(isoforms_bed)

    output:
        tuple val(sample), path("${chunk.baseName}.counts.txt"), path("${chunk.baseName}.map.txt"), emit: result

    script:
    // filter_transcriptome_align.py -t is the worker pool size; its parent process
    // needs one more core.
    def workers = Math.max(1, (task.cpus as int) - 1)
    """
    set -euo pipefail
    # filter_transcriptome_align.py is not installed in the env's bin/: it lives in
    # the flair package folder, which flair only adds to PATH while flair itself is
    # running. Call it by full path, the same way flair quantify ends up running it.
    FLAIR_PKG=\$(python3 -c "import os, flair; print(os.path.dirname(os.path.realpath(flair.__file__)))")
    test -f "\$FLAIR_PKG/filter_transcriptome_align.py"

    samtools view -h ${chunk} \\
    | python3 "\$FLAIR_PKG/filter_transcriptome_align.py" \\
        -s - \\
        -o ${chunk.baseName}.counts.txt \\
        -t ${workers} \\
        --quality 0 \\
        --trust_ends \\
        --check_splice \\
        -i ${isoforms_bed} \\
        --generate_map ${chunk.baseName}.map.txt

    # a chunk with no assigned reads still needs both files to exist
    touch ${chunk.baseName}.counts.txt ${chunk.baseName}.map.txt
    """

    stub:
    """
    echo -e "isoA_gene1\\t2\\nisoB_gene1\\t1" > ${chunk.baseName}.counts.txt
    echo -e "isoA_gene1\\t1:2|${chunk.baseName}_r1,1:2|${chunk.baseName}_r2\\nisoB_gene1\\t1:2|${chunk.baseName}_r3" > ${chunk.baseName}.map.txt
    """
}


process FLAIR_QUANT_MERGE {

    tag { sample }
    publishDir "${params.outdir}/preprocessing/flair/quantify/${sample}", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(counts, stageAs: 'counts/*'), path(maps, stageAs: 'maps/*')

    output:
        tuple val(sample), path("${sample}.flair.counts.tsv"),     emit: counts
        tuple val(sample), path("${sample}.isoform.read.map.txt"), emit: read_map

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
import glob, sys
from collections import Counter

sample = "${sample}"
count_files = sorted(glob.glob("counts/*"))
map_files   = sorted(glob.glob("maps/*"))
if not count_files or len(count_files) != len(map_files):
    sys.exit(f"[FLAIR_QUANT_MERGE] {sample}: {len(count_files)} counts vs {len(map_files)} map files")

# counts: isoform <tab> n, summed over chunks.
# Output matches flair quantify --sample_id_only: header 'ID<tab>sample', features sorted, integer counts.
counts = Counter()
for f in count_files:
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\\n").split("\\t")
            if len(p) >= 2 and p[0]:
                counts[p[0]] += int(float(p[1]))

with open(f"{sample}.flair.counts.tsv", "w") as out:
    out.write(f"ID\\t{sample}\\n")
    for iso in sorted(counts):
        out.write(f"{iso}\\t{counts[iso]}\\n")

# maps: concatenated (one row per isoform per chunk; build_read_isoform_lookup.py
# already accepts repeated isoform rows, as MERGE_FLAIR_MAPS produces them too).
# Sanity checks: map totals must equal counts; warn if any read is assigned twice.
seen = set()
dup = 0
n_map = 0
with open(f"{sample}.isoform.read.map.txt", "w") as out:
    for f in map_files:
        with open(f) as fh:
            for line in fh:
                if not line.strip():
                    continue
                out.write(line if line.endswith("\\n") else line + "\\n")
                p = line.rstrip("\\n").split("\\t")
                if len(p) < 2:
                    continue
                for r in p[1].split(","):
                    r = r.strip()
                    if not r:
                        continue
                    n_map += 1
                    if r in seen:
                        dup += 1
                    else:
                        seen.add(r)

n_counts = sum(counts.values())
print(f"[FLAIR_QUANT_MERGE] {sample}: {len(count_files)} chunks, {len(counts)} isoforms, "
      f"{n_counts} assigned reads, {n_map} map entries, {dup} duplicate read IDs")
if dup:
    # Each read is assigned at most once, so this means duplicate read names in the
    # input FASTQ. Warn rather than fail so a long run is not lost.
    print(f"[FLAIR_QUANT_MERGE] WARNING {sample}: {dup} read IDs assigned more than once "
          f"(duplicate read names in the input FASTQ?)", file=sys.stderr)
if n_map != n_counts:
    sys.exit(f"[FLAIR_QUANT_MERGE] {sample}: map entries ({n_map}) != summed counts ({n_counts})")
PY
    """
}
