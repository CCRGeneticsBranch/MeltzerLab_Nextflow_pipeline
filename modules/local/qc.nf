
process Fastqc {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/qc/", mode: 'copy',pattern: "fastqc"

    input:
    tuple val(meta), path(trim), path(r1fq), path(r2fq)

    output:
    tuple val(meta), path("fastqc") , emit: fastqc_results
    path "versions.yml"             , emit: versions


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    if [ ! -d fastqc ];then mkdir -p fastqc;fi
    fastqc --extract ${trim[0]} ${trim[1]} $r1fq $r2fq -t $task.cpus -o fastqc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Fastqc: \$(fastqc --version|awk '{print \$2}')
    END_VERSIONS
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

process RNAseQC {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/qc/${meta.lib}", mode: 'copy'

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(rRNA_interval),
        path(transcript_gtf)

    output:
    tuple val(meta),
        path("rnaseqc/report.html")

    stub:
     """
     touch "report.html"
     """

    script:
     """
     java -jar \$RNASEQCJAR -r ${genome} -rRNA ${rRNA_interval} -o rnaseqc -s "${meta.lib}|${bam}|${meta.lib}" -t ${transcript_gtf}

     """

}


process Kraken2 {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/qc/${meta.lib}/kraken", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(r1fq), path(r2fq),path(kraken2_db)

    output:
    tuple val(meta),path("${meta.lib}.kraken2_output.txt"), emit: kraken_output
    tuple val(meta),path("${meta.lib}_kraken2.report.txt"), emit : kraken_report
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.kraken_output"
    touch "${meta.lib}.report"

    """


    script:
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    kraken2 --db ${kraken2_db} --gzip-compressed --threads ${task.cpus} --output ${prefix}.kraken2_output.txt --paired ${r1fq} ${r2fq} --use-names --report ${prefix}_kraken2.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Kraken: \$(kraken2 --version|head -1|awk '{print \$3}')
    END_VERSIONS
    """
}
