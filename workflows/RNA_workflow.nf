
include {Fastp} from '../modules/local/Fastp.nf'
include {Star} from '../modules/local/star.nf'



workflow RNA_workflow {

star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
aligner                 = Channel.value(params.RNA_aligner)

take:
     samples_ch

main:

Fastp(samples_ch)

Star(Fastp.out
    .combine(star_genomeIndex)
    .combine(gtf)
    .combine(aligner)
)



}
