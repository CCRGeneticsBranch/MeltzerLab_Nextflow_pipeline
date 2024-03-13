process Picard_AddReadgroups {
        tag "$meta.lib"

        input:
        tuple val(meta), path(bam),path(index),val(aligner)

        output:
        tuple val(meta),
        path("${meta.lib}.${aligner}_${meta.genome}.bam"),
        path("${meta.lib}.${aligner}_${meta.genome}.bam.bai")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        set -exo pipefail

        java -Xmx10g  -jar \$PICARDJAR AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $bam  -OUTPUT ${prefix}.${aligner}_${meta.genome}.bam -SORT_ORDER coordinate -RGLB ${prefix} -RGPU ${prefix} -RGPL ILLUMINA -RGSM ${prefix} -RGCN khanlab

        java -Xmx10g  -jar \$PICARDJAR BuildBamIndex  -INPUT ${prefix}.${aligner}_${meta.genome}.bam -OUTPUT ${prefix}.${aligner}_${meta.genome}.bam.bai
        """
}

process Picard_MarkDuplicates {
        tag "$meta.lib"
        publishDir "${params.resultsdir}/qc/picard_metrics",mode: 'copy',pattern: "*markdup.txt"
        input:

        tuple val(meta),
        path(bam),
        path(index),
        val(aligner)

        output:
        tuple val(meta), path("${meta.lib}.${aligner}_${meta.genome}.dd.bam") , emit:bam
        tuple val(meta), path("${meta.lib}.${aligner}_${meta.genome}.dd.bam.bai") , emit:bai
        path("${meta.lib}.${aligner}_${meta.genome}.markdup.txt") , emit: markdup


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        java  -jar -Xmx10g \$PICARDJAR MarkDuplicates AS=true M=${prefix}.${aligner}_${meta.genome}.markdup.txt INPUT=$bam OUTPUT=${prefix}.${aligner}_${meta.genome}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar -Xmx10g \$PICARDJAR BuildBamIndex  -INPUT ${prefix}.${aligner}_${meta.genome}.dd.bam -OUTPUT ${prefix}.${aligner}_${meta.genome}.dd.bam.bai

        """
}
