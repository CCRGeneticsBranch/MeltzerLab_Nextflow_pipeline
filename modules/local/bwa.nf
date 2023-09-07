process BWA_mem2 {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/bam", mode: 'copy'
    
    input:
      tuple val(meta), path(trim), path(bwa_genomeindex), val(DNA_aligner)

    output:
      tuple val(meta), path("${meta.lib}.${meta.id}.${DNA_aligner}.${meta.genome}.bam")

    script:
    """
    
    bwa-mem2 mem  -t ${task.cpus} -M  -R '@RG\\tID:${meta.lib}.${meta.id}\\tSM:${meta.lib}.${meta.id}\\tLB:${meta.lib}.${meta.id}\\tPL:illumina' ${bwa_genomeindex}/${meta.genome} ${trim[0]} ${trim[1]} | samblaster | samtools view -Sb - > ${meta.lib}.${meta.id}.${DNA_aligner}.${meta.genome}.bam
          
    """

}