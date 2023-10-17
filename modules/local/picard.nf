process Picard_AddReadgroups {
        tag "$meta.lib"

        input:
        tuple val(meta), path(bam),path(index)

        output:
        tuple val(meta),
        path("${meta.lib}_star.bam"),
        path("${meta.lib}_star.bam.bai")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        set -exo pipefail

        java -Xmx10g  -jar \$PICARDJAR AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $bam  -OUTPUT ${prefix}_star.bam -SORT_ORDER coordinate -RGLB ${prefix} -RGPU ${prefix} -RGPL ILLUMINA -RGSM ${prefix} -RGCN khanlab

        java -Xmx10g  -jar \$PICARDJAR BuildBamIndex  -INPUT ${prefix}_star.bam -OUTPUT ${prefix}_star.bam.bai
        """
}

process Picard_MarkDuplicates {
        tag "$meta.lib"

        input:

        tuple val(meta),
        path(bam),
        path(index)

        output:
        tuple val(meta),
        path("${meta.lib}.dd.bam"),
        path("${meta.lib}.dd.bam.bai")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        java  -jar -Xmx10g \$PICARDJAR MarkDuplicates AS=true M=${prefix}.markdup.txt INPUT=$bam OUTPUT=${prefix}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar -Xmx10g \$PICARDJAR BuildBamIndex  -INPUT ${prefix}.dd.bam -OUTPUT ${prefix}.dd.bam.bai

        """
}
