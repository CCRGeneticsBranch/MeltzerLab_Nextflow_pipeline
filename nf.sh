#!/bin/bash

if [[ "$#" -ne "2" ]]; then
    echo "   "
    echo "Usage: ./nf.sh <profile> <samplesheet> [output_directory]"
    echo "   "
    echo "This script requires two positional arguments:"
    echo "1. <profile>: Specify the profile to use (e.g., hg19). More genome support coming soon."
    echo "2. <samplesheet>: Provide the path to the samplesheet in Meltzer lab format."
    echo "Optional Input:"
    echo "3. [output_directory]: Customize the output directory path (default: $WF_HOME/results)."
    echo "   "
    exit
fi

module load nextflow/23.04.3 singularity  graphviz

SCRIPT_NAME="$0"
SCRIPT_DIRNAME=$(readlink -f $(dirname $0))
SCRIPT_BASENAME=$(basename $0)
WF_HOME=$SCRIPT_DIRNAME

export PROFILE=$1
export SAMPLESHEET=$2

if [ -n "$3" ]; then
    OUTDIR="$3"
else
    OUTDIR="$WF_HOME/results"
fi

nf_cmd="nextflow"
nf_cmd="$nf_cmd run"
nf_cmd="$nf_cmd $WF_HOME/main.nf -resume "
nf_cmd="$nf_cmd -profile $PROFILE"
#nf_cmd="$nf_cmd -with-trace"
#nf_cmd="$nf_cmd -with-timeline"
nf_cmd="$nf_cmd --samplesheet $SAMPLESHEET"
nf_cmd="$nf_cmd --resultsdir $OUTDIR"

echo $nf_cmd

eval $nf_cmd
