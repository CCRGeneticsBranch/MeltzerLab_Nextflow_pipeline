process Fastp {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/qc/Fastp", mode: 'copy', pattern: ["*.html", "*.json"]

    input:
        tuple val(meta), path(r1fq), path(r2fq)
    output:
        tuple val(meta), path("Trimmed_Reads/*.gz") , emit: trim_reads
        path("${meta.lib}.${meta.id}-${meta.genome}.html") , emit: html
        path("${meta.lib}.${meta.id}-${meta.genome}.json") , emit: json


    script:
    """
    mkdir Trimmed_Reads
    fastp --in1 ${r1fq} \
          --in2 ${r2fq}  \
          --out1 Trimmed_Reads/trim_${r1fq} \
          --out2 Trimmed_Reads/trim_${r2fq} \
          --thread $task.cpus \
          --html ${meta.lib}.${meta.id}-${meta.genome}.html \
          --json ${meta.lib}.${meta.id}-${meta.genome}.json

    """
}
