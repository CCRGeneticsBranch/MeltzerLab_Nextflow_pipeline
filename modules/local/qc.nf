
process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

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

    input:
    tuple val(meta),path(bam),path(bai),val(aligner)

    output:
    tuple val(meta),path("${meta.lib}.${meta.id}.${aligner}.${meta.genome}.flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${meta.lib}.${meta.id}.${aligner}.${meta.genome}.flagstat.txt

    """

}

