#!/bin/bash

if [[ "$#" -ne 1 ]]; then
    echo "   "
    echo "Usage: ./nf.sh  <samplesheet> "
    echo "   "
    echo "This script requires one positional argument:"
    echo "<samplesheet>: Provide the path to the samplesheet in Meltzer lab format."
    exit
fi

module load nextflow/23.10.0 singularity  graphviz

SCRIPT_NAME="$0"
SCRIPT_DIRNAME=$(readlink -f $(dirname $0))
SCRIPT_BASENAME=$(basename $0)
WF_HOME=$SCRIPT_DIRNAME


CONFIG_FILE="$WF_HOME/nextflow.config"


export SAMPLESHEET=$1

export OUTDIR=$(dirname $SAMPLESHEET)

if [ -d "$OUTDIR/RESULTS" ]; then
    RESULTSDIR="$OUTDIR/RESULTS"
else
    RESULTSDIR="$OUTDIR/RESULTS"
    mkdir "$RESULTSDIR"
fi
export RESULTSDIR

export NXF_HOME="$OUTDIR/.nextflow"

printenv|grep NXF


export WORKDIR="$OUTDIR/work"
export INPUTDIR="$OUTDIR/fastq"

nf_cmd="nextflow"
nf_cmd="$nf_cmd run"
nf_cmd="$nf_cmd -c $CONFIG_FILE"
nf_cmd="$nf_cmd $WF_HOME/main.nf -resume "
nf_cmd="$nf_cmd -profile biowulf"
#nf_cmd="$nf_cmd -with-trace"
#nf_cmd="$nf_cmd -with-timeline"
nf_cmd="$nf_cmd --samplesheet $SAMPLESHEET"
nf_cmd="$nf_cmd --resultsdir $RESULTSDIR"
nf_cmd="$nf_cmd --input_dir $INPUTDIR"

echo $nf_cmd

eval $nf_cmd
