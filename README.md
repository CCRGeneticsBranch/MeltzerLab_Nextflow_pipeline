# Meltzerlab Data Processing Pipeline

This repository contains a Nextflow based DNA/RNA data processing pipeline built specifically for Meltzer lab. This pipeline is built to run variant filtering and extensive QC on a sequencing run. The pipeline supports processing data from multiple references. We currently support

| Genome | Version          |
|--------|------------------|
| Human  | hg19, hg38       |
| Mouse  | mm10, mm39       |
| Dog    | canFam3, canFam6 |

# Pipeline Overview

![Meltzarlab_dataProcessing_workflow (1)](https://github.com/CCRGeneticsBranch/MeltzerLab_Nextflow_pipeline/assets/8798661/0b60bb1e-4770-41f1-90be-368cffd71cc9)


This pipeline is developed and deployed solely on NIH Biowulf.

# Biowulf

This pipeline is hosted under /data/GBNCI directory. You'll need to start an interactive session inorder to launch the pipeline.

` sinteractive --mem=30g --cpus-per-task=4 `

## Usage

```
Usage: /data/GBNCI/MeltzerLab_Nextflow_pipeline/nf.sh  <samplesheet> 
   
This script requires one positional argument:
<samplesheet>: Provide the full path to the samplesheet in a Meltzer lab accepted format.
eg: /data/GBNCI/MeltzerLab_Nextflow_pipeline/nf.sh /data/GBNCI/DATA/VG_test/samplesheet.tsv 

```

## Input Samplesheet

Input file is a tab seperated samplesheet file from Meltzer lab. However the pipeline can be launched as long as these columns are present.
Here is an example of the expected format.

| flowcell | library       | sample_type | genome | read1                   | read2                    | capture_targets             |
|----------|---------------|-------------|--------|-------------------------|--------------------------|-----------------------------|
| AAZ123   | L080123_RNA   | Total RNA   | hg19   | L080123_RNA_R1.fastq.gz | L080123_RNA_R1.fastq.gz | APC_TP53_BAC_KlenowBait   |


## Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/TOOL_NAME/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/TOOL_NAME/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](docs/CONTRIBUTING.md).

## References

This repo was originally generated from the [CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).
The template takes inspiration from nektool[^1] and the nf-core template.
If you plan to contribute your pipeline to nf-core, don't use this template -- instead follow nf-core's instructions[^2].

[^1]: nektool https://github.com/beardymcjohnface/nektool
[^2]: instructions for nf-core pipelines https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core
