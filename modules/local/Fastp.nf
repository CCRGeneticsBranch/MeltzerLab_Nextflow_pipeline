process Fastp {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/qc/${meta.lib}.${meta.id}/Fastp", mode: 'copy'

    input:
        tuple val(meta), path(r1fq), path(r2fq)
    output:
        tuple val(meta), path("Trimmed_Reads/*.gz") , emit: trim_reads

    script:
    """
    mkdir Trimmed_Reads
    fastp --in1 ${r1fq} \
          --in2 ${r2fq}  \
          --out1 Trimmed_Reads/trim_${r1fq} \
          --out2 Trimmed_Reads/trim_${r2fq} \
          --thread $task.cpus \
          --html Trimmed_Reads/${meta.lib}.${meta.id}

    """
}
