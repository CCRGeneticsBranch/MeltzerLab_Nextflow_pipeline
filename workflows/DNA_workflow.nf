include {Fastp} from '../modules/local/Fastp.nf'
include {BWA_mem2} from '../modules/local/bwa.nf'

workflow DNA_workflow {

bwa_genomeindex = Channel.of(file(params.bwa_genomeindex, checkIfExists:true))
DNA_aligner     = Channel.value(params.DNA_aligner)


take: 
     samples_ch

main:

Fastp(samples_ch)

BWA_mem2(Fastp.out
         .combine(bwa_genomeindex)
         .combine(DNA_aligner)
)

}