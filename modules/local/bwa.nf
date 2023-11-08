process BWA_mem2 {
    tag "$meta.lib"
    //publishDir "${params.resultsdir}/bwa", mode: 'copy'

    input:
      tuple val(meta), path(trim), path(bwa_genomeindex), val(aligner)

    output:
      tuple val(meta),
      path("${meta.lib}.${meta.id}_${aligner}_${meta.genome}.bam"),
      path("${meta.lib}.${meta.id}_${aligner}_${meta.genome}.bam.bai")


    script:
    """

    bwa-mem2 mem  -t ${task.cpus} -M  -R '@RG\\tID:${meta.lib}.${meta.id}\\tSM:${meta.lib}.${meta.id}\\tLB:${meta.lib}.${meta.id}\\tPL:illumina' ${bwa_genomeindex}/${meta.genome} ${trim[0]} ${trim[1]} | samblaster -M | samtools view -Sb - |  samtools sort -m 30000000000 - > ${meta.lib}.${meta.id}_${aligner}_${meta.genome}.bam
    samtools index ${meta.lib}.${meta.id}_${aligner}_${meta.genome}.bam
    """

}
