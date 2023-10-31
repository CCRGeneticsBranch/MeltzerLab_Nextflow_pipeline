
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
        Kraken2
        Krona
        Multiqc} from '../modules/local/qc.nf'

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

Krona(Kraken2.out.kraken_output)


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
picard_output = Picard_MarkDuplicates.out.bam.combine(Picard_MarkDuplicates.out.bai,by:[0])
//picard_output.view()

GATK_SplitNCigarReads(picard_output
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
     picard_output
     .combine(genome)
     .combine(genome_fai)
     .combine(genome_dict)
     .combine(rRNA_interval)
     .combine(transcript_gtf)
)

multiqc_input = Fastqc.out.fastqc_results
               .join(Kraken2.out.kraken_report)
               .join(Krona.out.krona_output)
               .join(Flagstat.out)
               .join(Idxstats.out)
               .join(CollectMultipleMetrics.out)
               .join(Picard_MarkDuplicates.out.markdup)

multiqc_input_files = multiqc_input.map { tuple -> tuple.drop(1) }
multiqc_input_meta = multiqc_input.map { tuple -> tuple[0] }



Multiqc(multiqc_input_files,
           multiqc_input_meta)



}
