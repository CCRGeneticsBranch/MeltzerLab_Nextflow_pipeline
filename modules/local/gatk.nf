
process GATK_SplitNCigarReads {
        tag "$meta.lib"

        input:
        tuple val(meta),
        path(bam),
        path(index),
        path(ref_folder),
        val(aligner)

        output:
        tuple val(meta),
        path("${meta.lib}.${aligner}_${meta.genome}.split.bam"),
        path("${meta.lib}.${aligner}_${meta.genome}.split.bai")


        script:

        """
        java -Xmx70g -jar \$GATK_JAR SplitNCigarReads -R ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa -I $bam -O ${meta.lib}.${aligner}_${meta.genome}.split.bam

        """
}


process GATK_BaseRecalibrator {
    tag "$meta.lib"

    input:
      tuple val(meta),
      path(bam),
      path(bai),
      path(ref_folder),
      val(aligner)

    output:
      tuple val(meta),
      path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt")


    script:
    """
    if [[ "${meta.genome}" == "hg19" || "${meta.genome}" == "hg38" ]]; then
      java -Xmx70g -jar \$GATK_JAR  BaseRecalibrator \
        -R ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_dbsnp.vcf.gz \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_Mills_and_1000G_gold_standard.indels.vcf.gz \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_1000G_phase1.snps.high_confidence.vcf.gz \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_1000G_phase1.indels.vcf.gz \
        -I ${bam} \
        -O ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt
    elif [[ "${meta.genome}" == "mm10" || "${meta.genome}" == "mm39" ]]; then
      java -Xmx70g -jar \$GATK_JAR  BaseRecalibrator \
        -R ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_snp.vcf \
        --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_indels.vcf \
        -I ${bam} \
        -O ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt
    elif [[ "${meta.genome}" == "canFam3" || "${meta.genome}" == "canFam6" ]]; then
      java -Xmx70g -jar \$GATK_JAR  BaseRecalibrator \
      -R ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa \
      --known-sites  ${ref_folder}/${meta.genome}/${meta.genome}_snp.vcf \
      -I ${bam} \
      -O ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt
    fi
    """

}

process GATK_ApplyBQSR {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/bam", mode: 'copy'

    input:
      tuple val(meta),
      path(bam),
      path(bai),
      path(recalibration),
      path(ref_folder),
      val(aligner)

    output:
      tuple val(meta),
      path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.final.bam"),
      path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.final.bai")


    script:
    """
    java -Xmx70g -jar \$GATK_JAR  ApplyBQSR \
      -R ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa \
      -I ${bam} \
      --bqsr-recal-file ${recalibration} \
      -O ${meta.lib}.${meta.id}.${aligner}-${meta.genome}.final.bam

    """

}
