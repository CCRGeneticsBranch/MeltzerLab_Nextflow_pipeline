
process Fastqc {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/qc/fastqc", mode: 'copy',pattern: "fastqc"

    input:
    tuple val(meta), path(trim), path(r1fq), path(r2fq)

    output:
    tuple val(meta), path("${meta.lib}_fastqc") , emit: fastqc_results
    path "versions.yml"             , emit: versions


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    if [ ! -d ${meta.lib}_fastqc ];then mkdir -p ${meta.lib}_fastqc;fi
    fastqc --extract ${trim[0]} ${trim[1]} $r1fq $r2fq -t $task.cpus -o ${meta.lib}_fastqc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Fastqc: \$(fastqc --version|awk '{print \$2}')
    END_VERSIONS
    """
}

process Flagstat {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/qc/samtools", mode: 'copy'
    input:
    tuple val(meta),path(bam),path(bai),val(aligner)

    output:
    tuple val(meta),path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.flagstat.tsv")

    script:
    """
    samtools flagstat --output-fmt tsv ${bam} > ${meta.lib}.${meta.id}.${aligner}-${meta.genome}.flagstat.tsv

    """

}

process Idxstats {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/qc/samtools", mode: 'copy'

    input:
    tuple val(meta),path(bam),path(bai),val(aligner)

    output:
    tuple val(meta),path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.idxstats.tsv")

    script:
    """
    samtools idxstats ${bam} > ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.idxstats.tsv

    """

}

process CollectMultipleMetrics {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/qc/picard_metrics", mode: 'copy'
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
    publishDir "${params.resultsdir}/qc/${meta.lib}", mode: 'copy'

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

    publishDir "${params.resultsdir}/qc/kraken", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(r1fq), path(r2fq),path(kraken2_db)

    output:
    tuple val(meta),path("${meta.lib}-${meta.genome}.kraken2_output.txt"), emit: kraken_output
    tuple val(meta),path("${meta.lib}-${meta.genome}.kraken2.report.txt"), emit : kraken_report
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.kraken_output"
    touch "${meta.lib}.report"

    """


    script:
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    kraken2 --db ${kraken2_db} --gzip-compressed --threads ${task.cpus} --output ${prefix}-${meta.genome}.kraken2_output.txt --paired ${r1fq} ${r2fq} --use-names --report ${prefix}-${meta.genome}.kraken2.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Kraken: \$(kraken2 --version|head -1|awk '{print \$3}')
    END_VERSIONS
    """
}


process Krona {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/qc/kraken", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(kraken2_output)

    output:
    tuple val(meta),path("${meta.lib}-${meta.genome}.kraken2.krona.html"), emit: krona_output


    stub:
    """
    touch "${meta.lib}-${meta.genome}.kraken2.krona.html"

    """


    script:
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    ktImportTaxonomy -q 2 -t 3 ${kraken2_output} -o ${prefix}-${meta.genome}.kraken2.krona.html

    """
}


process Multiqc {
    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.lib}/qc", mode: 'copy',pattern: "*html"

    input:
    path(input_files)
    val(meta)


    output:
    path("multiqc_report.html") , emit: multiqc_report
    path "versions.yml"             , emit: versions

    script:
    """

    echo  "${input_files.join('\n')}" > multiqc_input_files
    multiqc --file-list multiqc_input_files -f

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    multiqc: \$(multiqc --version|sed 's/.*version //')
END_VERSIONS
    """

}
