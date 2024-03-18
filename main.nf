nextflow.enable.dsl=2

log.info """\
TOOL_NAME
=============
NF version   : $nextflow.version
runName      : $workflow.runName
username     : $workflow.userName
configs      : $workflow.configFiles
profile      : $workflow.profile
cmd line     : $workflow.commandLine
start time   : $workflow.start
projectDir   : $workflow.projectDir
launchDir    : $workflow.launchDir
workDir      : $workflow.workDir
homeDir      : $workflow.homeDir
"""
.stripIndent()


include {Check_Input} from './modules/local/check_input.nf'
include {DNA_workflow} from './workflows/DNA_workflow'
include {RNA_workflow} from './workflows/RNA_workflow'
include {NGSCheckMate
        Ncm_data_processing
        Multiqc} from './modules/local/qc.nf'


workflow  {


Check_Input(file(params.samplesheet, checkIfExists: true))


 def input_dir = params.input_dir
 data = Check_Input.out
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        def meta = [:]
        meta.id = row.flowcell
        meta.lib = row.library
        meta.sc = row.capture_targets
        meta.type = row.sample_type
        meta.genome = row.genome
        def read1Path = "${input_dir}/${row.read1}"
        def read2Path = "${input_dir}/${row.read2}"

        def fastq_meta = [meta,read1Path, read2Path]

        return fastq_meta
    }

samples = data.branch {
    RNA: it[0].type.contains("RNA")
    DNA: it[0].type.contains("DNA")
}


samples.RNA | RNA_workflow
samples.DNA | DNA_workflow

vaf_ch = DNA_workflow.out.ncm_vaf.ifEmpty([]).merge(RNA_workflow.out.ncm_vaf).ifEmpty([])
        .filter { it.size() > 0 }

vaf_ch|NGSCheckMate

NGSCheckMate.out.pdf|Ncm_data_processing


multiqc_ch = DNA_workflow.out.multiqc_ch.collect().ifEmpty([]).merge(RNA_workflow.out.multiqc_ch.collect().ifEmpty([]))
            .filter { it.size() > 0 }

multiqc_ch|Multiqc

}
