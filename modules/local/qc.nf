
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

process Fastq_screen {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/qc/fastq_screen", mode: 'copy'

    input:
    tuple val(meta),path(trim),path(fastq_screen_config),path(fqs_db)

    stub:
    """
    touch "${meta.lib}_R1_screen.html"
    touch "${meta.lib}_R2_screen.html"
    """

    output:
    tuple val(meta),path("*html"),path("*png"),path("*txt")


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    if [ ! -d fastq_screen ];then mkdir -p fastq_screen;fi
    ls ${fqs_db}
    fastq_screen --conf ${fastq_screen_config} --subset 1000000 --aligner bowtie2 --force ${trim[0]} ${trim[1]}
    """
}

process NGSCheckMate_vaf {
    tag "$meta.lib"
    // Running two step ngscheckmate. At library level we generate vaf file. later we run step 2 at run level.
    publishDir "${params.resultsdir}/qc/ncm/vaf", mode: 'copy'

    input:
    tuple val(meta),path(trim),val(aligner)

    stub:
    """
    touch "${meta.lib}.${meta.id}.${aligner}-${meta.genome}.vaf"
    """

    output:
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.vaf")


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT

    \$NCM_HOME/ngscheckmate_fastq -p ${task.cpus} \
            -1 ${trim[0]} \
            -2 ${trim[1]} \
            \$NCM_HOME/SNP/SNP.pt \
            > ${meta.lib}.${meta.id}.${aligner}-${meta.genome}.vaf
    """
}


process NGSCheckMate {

    publishDir "${params.resultsdir}/qc/ncm", mode: 'copy'

    input:
    path(vaf_files)

    stub:
    """
    touch "NGSCheckMate.pdf"
    touch "NGSCheckMate_all.txt"
    """

    output:
    path("NGSCheckMate.pdf") , emit : pdf
    path("NGSCheckMate_all.txt"), emit : png


    script:
    def args = task.ext.args   ?: ''

    """
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT
    mv *vaf \$TMP

    python \$NCM_HOME/vaf_ncm.py -I \$TMP \
            -O \$TMP \
            -N NGSCheckMate
    cp \$TMP/NGSCheckMate.pdf .
    echo -e "Sample1\tmatched/unmatched\tSample2\tCorrelation\tDepth"| cat - \$TMP/NGSCheckMate_all.txt > NGSCheckMate_all.txt

    """
}

process Ncm_data_processing{

    publishDir "${params.resultsdir}/qc/ncm", mode: 'copy'

    input:
    path(ncm_pdf)

    stub:
    """
    touch "NGSCheckMate.png"
    """

    output:
    path("NGSCheckMate.png")

    script:
    def args = task.ext.args   ?: ''

    """
    pdftoppm -png -r 600 ${ncm_pdf} -singlefile NGSCheckMate

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
    path(ref_folder),
    val(aligner)

    output:
    tuple val(meta),
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.quality_distribution_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.alignment_summary_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.insert_size_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.gc_bias.summary_metrics"),
    path("${meta.lib}.${meta.id}.${aligner}-${meta.genome}.quality_yield_metrics")

    script:
    """
    java -Xmx60g -jar \$PICARDJAR CollectMultipleMetrics VALIDATION_STRINGENCY=SILENT \
    INPUT=${bam} \
    OUTPUT=${meta.lib}.${meta.id}.${aligner}-${meta.genome} \
    REFERENCE_SEQUENCE=${ref_folder}/${meta.genome}/Index_files/${meta.genome}.fa \
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
