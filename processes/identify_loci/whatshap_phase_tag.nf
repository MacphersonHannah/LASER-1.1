// Takes the output vcf from DeepVariant and the haplotype-resolved bam from Whatshap to produce a phased VCF and haplotagged BAM.
// to see what it's actually used to phase, also runs whatshap stats phased.vcf.gz.
// 01/06/2026 added handling for samples with no PASS heterozygous variants after filtering (eg low-coverage barcodes). These produce an empty phased VCF and an untagged BAM
// rather than crashing the pipeline. Prints a log warning if so
process WHATSHAP_PHASE_TAG {
    tag { sample }
    publishDir "${params.outdir}/identify_loci/whatshap/${sample}", mode: 'link'
    conda "${projectDir}/envs_prebuilt/core_utils_env"

    input:
    tuple val(sample), path(vcf), path(bam), path(ref_fa), path(ref_fai)

    output:
    tuple val(sample), path("phased.vcf.gz"), emit: phased_vcf
    tuple val(sample), path("${sample}.flnc.hp.bam"), emit: hp_bam
    path "${sample}.flnc.hp.bam.bai"
    path "${sample}.whatshap.stats.txt"
    path "phased.vcf.gz.tbi"

    script:
    def regionArg = params.hp_regions ? "--regions-file ${params.hp_regions}" : ""
    def indelFlag = params.phase_indels ? "--indels" : ""
    def minQual = params.phase_min_qual ?: 0
    def qualExpr = "QUAL>=${minQual}"

    """
    set -euo pipefail
    export LC_ALL=C

    # Ensure indexes exist
    if [ ! -s "${ref_fai}" ]; then samtools faidx "${ref_fa}"; fi
    if [ ! -s "${bam}.bai" ]; then samtools index "${bam}"; fi
    if [ ! -s "${vcf}.tbi" ]; then tabix -p vcf "${vcf}"; fi

    # Filter for WhatsHap: keep PASS heterozygous biallelic variants
    bcftools view \\
        -f PASS \\
        -g het \\
        -m2 -M2 \\
        -i "${qualExpr}" \\
        -Oz -o sites.biallelic.vcf.gz ${vcf}

    tabix -p vcf sites.biallelic.vcf.gz

    # Check if there are any variants to phase
    N_VARS=\$(bcftools view --no-header sites.biallelic.vcf.gz 2>/dev/null | wc -l || echo 0)

    if [ "\$N_VARS" -eq 0 ]; then
        echo "WARNING: ${sample} has no PASS heterozygous variants after filtering." >&2
        echo "WARNING: This is likely due to low coverage or no callable heterozygous SNVs." >&2
        echo "WARNING: Skipping phasing for ${sample} - output BAM will be untagged." >&2

        # Create empty phased VCF with valid header
        bcftools view -h ${vcf} | bgzip > phased.vcf.gz
        tabix -p vcf phased.vcf.gz

        # Copy original BAM as untagged output
        cp ${bam} ${sample}.flnc.hp.bam
        samtools index -@ ${task.cpus} ${sample}.flnc.hp.bam

        # Empty stats file
        echo "No variants to phase for ${sample}" > ${sample}.whatshap.stats.txt

    else
        # 3) PHASE
        whatshap phase \\
            --reference "${ref_fa}" \\
            --ignore-read-groups \\
            --distrust-genotypes \\
            ${indelFlag} \\
            -o phased.vcf.gz \\
            sites.biallelic.vcf.gz \\
            "${bam}"

        tabix -p vcf phased.vcf.gz

        # 4) HAPLOTAG
        whatshap haplotag \\
            --no-reference \\
            --ignore-read-groups \\
            --tag-supplementary \\
            --skip-missing-contigs \\
            ${regionArg} \\
            -o "${sample}.flnc.hp.bam" \\
            phased.vcf.gz \\
            "${bam}"

        # 5) STATS
        whatshap stats phased.vcf.gz > "${sample}.whatshap.stats.txt"
        samtools index -@ ${task.cpus} "${sample}.flnc.hp.bam"
    fi
    """
}
