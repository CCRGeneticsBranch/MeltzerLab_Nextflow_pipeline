include {Fastp} from '../modules/local/Fastp.nf'

workflow DNA_workflow {

take: 
     samples_ch

main:

Fastp(samples_ch)


}