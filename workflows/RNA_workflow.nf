
include {Fastp} from '../modules/local/Fastp.nf'
include {Star} from '../modules/local/star.nf'
include {Picard_AddReadgroups
        Picard_MarkDuplicates} from '../modules/local/picard.nf'
include {GATK_SplitNCigarReads
        GATK_BaseRecalibrator
        GATK_ApplyBQSR} from '../modules/local/gatk.nf'

include {Fastq_screen
        NGSCheckMate_vaf
        Flagstat
        Idxstats
        CollectMultipleMetrics
        Fastqc
        Kraken2
        Krona
        Multiqc
        RNAseQC
        Strandedness
        CollectRnaSeqMetrics} from '../modules/local/qc.nf'

workflow RNA_workflow {

kraken2_db = Channel.of(file(params.kraken2_db, checkIfExists:true))
ref_folder = Channel.of(file(params.ref_folder, checkIfExists:true))
fastq_screen_config = Channel.of(file(params.fastq_screen_config, checkIfExists:true))
fastq_screen_db = Channel.of(file(params.fastq_screen_db, checkIfExists:true))
RNA_aligner                 = Channel.value(params.RNA_aligner)

/*
star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
aligner                 = Channel.value(params.RNA_aligner)
genome                  = Channel.of(file(params.genome, checkIfExists:true))
genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
phase1_1000g            = Channel.of(file(params.phase1_1000g, checkIfExists:true))
Mills_and_1000g         = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))
phase1_1000g_idx            = Channel.of(file(params.phase1_1000g_idx, checkIfExists:true))
Mills_and_1000g_idx         = Channel.of(file(params.Mills_and_1000g_idx, checkIfExists:true))
rRNA_interval               = Channel.of(file(params.rRNA_interval, checkIfExists:true))
transcript_gtf              = Channel.of(file(params.transcript_gtf, checkIfExists:true))
*/
take:
     samples_ch

main:

multiqc_ch = Channel.of()
Fastqc(samples_ch)

multiqc_ch = multiqc_ch.mix(Fastqc.out.fastqc_zip)
multiqc_ch = multiqc_ch.mix(Fastqc.out.fastqc_html)
Kraken2(samples_ch
     .combine(kraken2_db))

Krona(Kraken2.out.kraken_output)


Fastp(samples_ch)
multiqc_ch = multiqc_ch.mix(Fastp.out.html)
multiqc_ch = multiqc_ch.mix(Fastp.out.json)
Fastq_screen_input = samples_ch
                        .combine(fastq_screen_config)
                        .combine(fastq_screen_db)

Fastq_screen(Fastq_screen_input)
multiqc_ch =multiqc_ch.mix(Fastq_screen.out)

check_genome = Fastp.out.trim_reads.branch {
     human: it[0].genome == "hg19" || it[0].genome == "hg38"
}

check_genome.human.combine(RNA_aligner)| NGSCheckMate_vaf


Star(Fastp.out.trim_reads
    .combine(ref_folder)
    .combine(RNA_aligner)
)

Strandedness(Star.out.genome_bam
     .join(Star.out.genome_bai,by:[0])
     .combine(ref_folder)
     .combine(RNA_aligner))

CollectRnaSeqMetrics(
     Star.out.genome_bam
     .join(Star.out.genome_bai,by:[0])
     .join(Strandedness.out,by:[0])
     .combine(ref_folder)
     .combine(RNA_aligner))

multiqc_ch =multiqc_ch.mix(CollectRnaSeqMetrics.out)

Picard_AddReadgroups(Star.out.genome_bam
          .join(Star.out.genome_bai,by:[0])
          .combine(RNA_aligner))

Picard_MarkDuplicates(Picard_AddReadgroups.out
          .combine(RNA_aligner))

multiqc_ch =multiqc_ch.mix(Picard_MarkDuplicates.out.markdup)

picard_output = Picard_MarkDuplicates.out.bam.combine(Picard_MarkDuplicates.out.bai,by:[0])


GATK_SplitNCigarReads(picard_output
          .combine(ref_folder)
          .combine(RNA_aligner))

GATK_BaseRecalibrator(
     GATK_SplitNCigarReads.out
     .combine(ref_folder)
     .combine(RNA_aligner)
)

Applybqsr_input = GATK_SplitNCigarReads.out.join(GATK_BaseRecalibrator.out,by:[0])

GATK_ApplyBQSR(
     Applybqsr_input
     .combine(ref_folder)
     .combine(RNA_aligner)
)

Flagstat(
     GATK_ApplyBQSR.out
     .combine(RNA_aligner)
)
multiqc_ch =multiqc_ch.mix(Flagstat.out)
Idxstats(
     GATK_ApplyBQSR.out
     .combine(RNA_aligner)
)
multiqc_ch =multiqc_ch.mix(Idxstats.out)
CollectMultipleMetrics(
     GATK_ApplyBQSR.out
     .combine(ref_folder)
     .combine(RNA_aligner)
)
multiqc_ch =multiqc_ch.mix(CollectMultipleMetrics.out)
RNAseQC(
     GATK_ApplyBQSR.out
     .combine(ref_folder)
     .combine(RNA_aligner)
)
multiqc_ch =multiqc_ch.mix(RNAseQC.out)

emit:
ncm_vaf = NGSCheckMate_vaf.out.collect()
multiqc_ch = multiqc_ch
}
