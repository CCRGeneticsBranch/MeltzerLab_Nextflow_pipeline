
process GATK_SplitNCigarReads {
        tag "$meta.lib"

        input:
        tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

        output:
        tuple val(meta),
        path("${meta.lib}.split.bam"),
        path("${meta.lib}.split.bai")


        script:

        """
        java -Xmx70g -jar \$GATK_JAR SplitNCigarReads -R $genome -I $bam -O ${meta.lib}.split.bam

        """
}


process GATK_BaseRecalibrator {
    tag "$meta.lib"

    input:
      tuple val(meta),
      path(bam),
      path(bai),
      path(genome),
      path(genome_fai),
      path(genome_dict),
      path(phase1_1000g),
      path(phase1_1000g_idx),
      path(Mills_and_1000g),
      path(Mills_and_1000g_idx),
      val(aligner)

    output:
      tuple val(meta),
      path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt")


    script:
    """
    java -Xmx70g -jar \$GATK_JAR  BaseRecalibrator \
      -R ${genome} \
      --known-sites  ${Mills_and_1000g} \
      --known-sites  ${phase1_1000g} \
      -I ${bam} \
      -O ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.recalibration.matrix.txt

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
      path(genome),
      path(genome_fai),
      path(genome_dict),
      val(aligner)

    output:
      tuple val(meta),
      path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}_final.bam"),
      path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}_final.bai")


    script:
    """
    java -Xmx70g -jar \$GATK_JAR  ApplyBQSR \
      -R ${genome} \
      -I ${bam} \
      --bqsr-recal-file ${recalibration} \
      -O ${meta.lib}.${meta.id}.${aligner}.${meta.genome}_final.bam

    """

}
