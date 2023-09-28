
process FASTQC {
    tag { sample_id }
    publishDir "${params.resultsdir}/${meta.id}/qc", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}*.html")

    script:
    """
    fastqc \
        $fastq \
        -t $task.cpus \
        -o .
    """
}

process Flagstat {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/qc", mode: 'copy'
    input:
    tuple val(meta),path(bam),path(bai),val(aligner)

    output:
    tuple val(meta),path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.flagstat.txt

    """

}

process Idxstats {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/qc", mode: 'copy'
    
    input:
    tuple val(meta),path(bam),path(bai),val(aligner)

    output:
    tuple val(meta),path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.idxstats.txt")

    script:
    """
    samtools idxstats ${bam} > ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.idxstats.txt

    """

}

process CollectMultipleMetrics {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/qc/picard_metrics", mode: 'copy'
    input:
    tuple val(meta),
    path(bam),
    path(bai),
    path(genome),
    path(genome_fai),
    path(genome_dict),
    val(aligner)

    output:
    tuple val(meta),
    path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.quality_distribution_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.alignment_summary_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.insert_size_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.gc_bias.summary_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.quality_yield_metrics")

    script:
    """
    java -Xmx60g -jar \$PICARDJAR CollectMultipleMetrics VALIDATION_STRINGENCY=SILENT \
    INPUT=${bam} \
    OUTPUT=${meta.lib}.${meta.id}.${aligner}.${meta.genome} \
    REFERENCE_SEQUENCE=${genome} \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=QualityScoreDistribution \
    PROGRAM=MeanQualityByCycle \
    PROGRAM=CollectBaseDistributionByCycle \
    PROGRAM=CollectGcBiasMetrics  \
    PROGRAM=CollectSequencingArtifactMetrics \
    PROGRAM=CollectQualityYieldMetrics  

    """


}