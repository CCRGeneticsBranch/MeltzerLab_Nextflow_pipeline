process Star {
    tag "$meta.lib"


    input:
    tuple val(meta), path(trim),path(ref_folder),val(aligner)

    output:
    tuple val(meta), path("${meta.lib}.Aligned.toTranscriptome.out.${aligner}_${meta.genome}.bam") , emit: transcriptome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam"), emit: genome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam.bai"), emit: genome_bai
    tuple val(meta), path("${meta.lib}.Chimeric.out.junction"), emit: chimeric_junction
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.Aligned.toTranscriptome.out.${aligner}_${meta.genome}.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam.bai"
    touch "${meta.lib}.Chimeric.out.junction"
    """
    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT


        # run STAR alignment
        STAR --genomeDir ${ref_folder}/${meta.genome}/Index_files/STAR_index/ \
            --readFilesIn ${trim[0]} ${trim[1]} \
            --readFilesCommand zcat \
            --sjdbGTFfile ${ref_folder}/${meta.genome}/Index_files/${meta.genome}.gtf \
            --runThreadN ${task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix ${prefix}. \
            --chimSegmentMin 12 \
            --chimOutJunctionFormat 1 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --chimSegmentReadGapMax 3 \
            --outFilterMismatchNmax 2 \
            --outSAMtype BAM Unsorted \
            --quantMode TranscriptomeSAM

        # sort files
        samtools sort -@ ${task.cpus}  -T \$TMP -o ${prefix}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam -O BAM ${prefix}.Aligned.out.bam


    # index files
    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.${aligner}_${meta.genome}.bam

    mv ${meta.lib}.Aligned.toTranscriptome.out.bam ${meta.lib}.Aligned.toTranscriptome.out.${aligner}_${meta.genome}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version)
    END_VERSIONS
    """
}
