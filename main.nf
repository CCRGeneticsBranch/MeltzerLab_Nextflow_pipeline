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



include {DNA_workflow} from './workflows/DNA_workflow'
include {RNA_workflow} from './workflows/RNA_workflow'
workflow  {


 def input_dir = params.input_dir
 data = Channel.of(file(params.samplesheet, checkIfExists:true))
 //data = Channel.fromPath("/data/GBNCI/DATA/AAC37LMM5/samplesheet.tsv")
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        def meta = [:]
        meta.id = row.flowcell
        meta.lib = row.library
        meta.sc = row.capture_targets
        meta.genome = row.genome
        meta.type = row.sample_type
        def read1Path = "${input_dir}/${row.read1}"
        def read2Path = "${input_dir}/${row.read2}"

        def fastq_meta = [meta,read1Path, read2Path]

        return fastq_meta

    }

samples = data.branch {
RNA: it[0].type == "Total RNA"
Tumor: it[0].type == "ChIP DNA"
}

samples.RNA | RNA_workflow
samples.Tumor | DNA_workflow


}
