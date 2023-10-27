
include {Fastp} from '../modules/local/Fastp.nf'
include {Star} from '../modules/local/star.nf'
include {Picard_AddReadgroups
        Picard_MarkDuplicates} from '../modules/local/picard.nf'
include {GATK_SplitNCigarReads
        GATK_BaseRecalibrator
        GATK_ApplyBQSR} from '../modules/local/gatk.nf'

include {Flagstat
        Idxstats
        CollectMultipleMetrics
        Fastqc
        RNAseQC
        Kraken2} from '../modules/local/qc.nf'

workflow RNA_workflow {

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

take:
     samples_ch

main:

kraken2_db = Channel.of(file(params.kraken2_db, checkIfExists:true))

Kraken2(samples_ch
     .combine(kraken2_db))


Fastp(samples_ch)

fastqc_input = Fastp.out.trim_reads.join(samples_ch, by:[0])

Fastqc(fastqc_input)

Star(Fastp.out
    .combine(star_genomeIndex)
    .combine(gtf)
    .combine(aligner)
)

Picard_AddReadgroups(Star.out.genome_bam.combine(Star.out.genome_bai,by:[0]))
Picard_MarkDuplicates(Picard_AddReadgroups.out)

GATK_SplitNCigarReads(Picard_MarkDuplicates.out
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
)

GATK_BaseRecalibrator(
     GATK_SplitNCigarReads.out
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
     .combine(phase1_1000g)
     .combine(phase1_1000g_idx)
     .combine(Mills_and_1000g)
     .combine(Mills_and_1000g_idx)
     .combine(aligner)
)


Applybqsr_input = GATK_SplitNCigarReads.out.join(GATK_BaseRecalibrator.out,by:[0])

GATK_ApplyBQSR(
     Applybqsr_input
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
     .combine(aligner)
)

Flagstat(
     GATK_ApplyBQSR.out
     .combine(aligner)
)

Idxstats(
     GATK_ApplyBQSR.out
     .combine(aligner)
)

CollectMultipleMetrics(
     GATK_ApplyBQSR.out
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
     .combine(aligner)
)
RNAseQC(
     Picard_MarkDuplicates.out
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
     .combine(rRNA_interval)
     .combine(transcript_gtf)
)

}
