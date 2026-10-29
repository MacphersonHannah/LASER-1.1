process RUN_CLAIR3_RNA {
    tag "${sample}"
    publishDir "${params.outdir}/preprocessing/clair3/${sample}", mode: 'link'

    container "${projectDir}/containers/clair3_rna_latest.sif"

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)
    val   model_name

    output:
    tuple val(sample), path("${sample}.raw.vcf.gz"), emit: vcf

    script:
    def phasingFlag = params.clair3_enable_phasing ? '--enable_phasing_model' : ''
    def minCov      = params.clair3_min_coverage != null ? params.clair3_min_coverage : 2
    """
    # setting temp directory within the work_ directory
    export TMPDIR=\$PWD

    /opt/bin/run_clair3_rna \\
        -B ${bam}                        \\
        -R ${ref_fa}                     \\
        -o clair3_rna_out                \\
        -t ${task.cpus}                  \\
        -p ${model_name}                 \\
        -P /opt/conda/envs/clair3_rna/bin/clair3_rna_models/${model_name}/variables \\
        --min_coverage ${minCov}         \\
        ${phasingFlag}

    # Fail if Clair3 silently produces nothing. Doesn't necessarily mean pipeline failed as it may genuinely have found no variants if very low read counts or a cell line or something, but means you can't do get any ASE analysis
    n_variants=\$(zgrep -c -v "^#" clair3_rna_out/output.vcf.gz || true)
    if [ "\$n_variants" -eq 0 ]; then
        echo "ERROR: Clair3-RNA called zero variants. Check container TMPDIR visibility and BAM content." >&2
        exit 1
    fi

    cp clair3_rna_out/output.vcf.gz ${sample}.raw.vcf.gz
    """

    stub:
    """
    touch ${sample}.raw.vcf.gz
    """
}