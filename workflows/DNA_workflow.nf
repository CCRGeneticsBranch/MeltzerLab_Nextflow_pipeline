include {Fastp} from '../modules/local/Fastp.nf'
include {BWA_mem2} from '../modules/local/bwa.nf'
include {GATK_BaseRecalibrator} from '../modules/local/gatk.nf'
include {GATK_ApplyBQSR} from '../modules/local/gatk.nf'
include {Fastq_screen
        NGSCheckMate_vaf
        Flagstat
        Idxstats
        CollectMultipleMetrics
        Fastqc
        Kraken2
        Krona
        Multiqc} from '../modules/local/qc.nf'

workflow DNA_workflow {

kraken2_db = Channel.of(file(params.kraken2_db, checkIfExists:true))
fastq_screen_config = Channel.of(file(params.fastq_screen_config, checkIfExists:true))
fastq_screen_db = Channel.of(file(params.fastq_screen_db, checkIfExists:true))
ref_folder = Channel.of(file(params.ref_folder, checkIfExists:true))
aligner     = Channel.value(params.DNA_aligner)

take:
     samples_ch

main:



Kraken2(samples_ch
     .combine(kraken2_db))

Krona(Kraken2.out.kraken_output)

Fastp(samples_ch)

fastqc_input = Fastp.out.trim_reads.join(samples_ch, by:[0])

Fastqc(fastqc_input)

Fastq_screen_input = Fastp.out.trim_reads
                        .combine(fastq_screen_config)
                        .combine(fastq_screen_db)

Fastq_screen(Fastq_screen_input)

check_genome = Fastp.out.trim_reads.branch {
     human: it[0].genome == "hg19" || it[0].genome == "hg38"
}

check_genome.human.combine(aligner)| NGSCheckMate_vaf



BWA_mem2(Fastp.out.trim_reads
         .combine(ref_folder)
         .combine(aligner)
)

GATK_BaseRecalibrator(
     BWA_mem2.out
     .combine(ref_folder)
     .combine(aligner)
)
Applybqsr_input = BWA_mem2.out.join(GATK_BaseRecalibrator.out,by:[0])
GATK_ApplyBQSR(
     Applybqsr_input
     .combine(ref_folder)
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
     .combine(ref_folder)
     .combine(aligner)
)
/*
multiqc_input = Fastqc.out.fastqc_results
               .join(Kraken2.out.kraken_report)
               .join(Krona.out.krona_output)
               .join(Flagstat.out)
               .join(Idxstats.out)
               .join(CollectMultipleMetrics.out)

multiqc_input_files = multiqc_input.map { tuple -> tuple.drop(1) }
multiqc_input_meta = multiqc_input.map { tuple -> tuple[0] }

Multiqc(multiqc_input_files,
           multiqc_input_meta)
*/
emit:
ncm_vaf = NGSCheckMate_vaf.out.collect()
}
