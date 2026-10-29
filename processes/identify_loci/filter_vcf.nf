// Filter a raw VCF to PASS variants only and index with tabix
// Runs outside the Clair3-RNA container using the whatshap conda env, which includes bcftools.
// Handles the case where Clair3 produces a valid but empty VCF (eg barcode with very low coverage / no callable variants). bcftools view fails on empty gzipped VCFs, so check first and create a proper empty output if needed.
nextflow.enable.dsl = 2

process FILTER_VCF {
    tag "${sample}"
    publishDir "${params.outdir}/preprocessing/clair3/${sample}", mode: 'link'
    conda "${projectDir}/envs_prebuilt/core_utils_env"

    input:
    tuple val(sample), path(raw_vcf)

    output:
    tuple val(sample), path("${sample}.clair3_rna.vcf.gz"),
                       path("${sample}.clair3_rna.vcf.gz.tbi"), emit: vcf

    script:
    """
    # Check if the VCF has any variant records (lines not starting with #)
    N_VARS=\$(bcftools view --no-header ${raw_vcf} 2>/dev/null | wc -l || echo 0)

    if [ "\$N_VARS" -eq 0 ]; then
        echo "WARNING: ${sample} VCF contains no variants - creating empty filtered VCF." >&2
        printf '##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' | \\
            bgzip > ${sample}.clair3_rna.vcf.gz
    else
        bcftools view \\
            --apply-filters PASS \\
            --output-type z \\
            --output ${sample}.clair3_rna.vcf.gz \\
            ${raw_vcf}
    fi

    bcftools index --tbi ${sample}.clair3_rna.vcf.gz
    """

    stub:
    """
    touch ${sample}.clair3_rna.vcf.gz
    touch ${sample}.clair3_rna.vcf.gz.tbi
    """
}
