<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Set up environment](#set-up-environment)
* [Run pipeline](#run-pipeline)

<!-- vim-markdown-toc -->

## Description

This is a Snakemake pipeline for RNAseq analysis with the primary aim of
comparing two groups of samples. Sample information is in `sample_sheet.tsv`.

## Set up environment

Assuming you have conda installed and configured for
[bioconda](https://bioconda.github.io/user/install.html) and you also have
[mamba](https://github.com/mamba-org/mamba) also installed. Create a dedicated
environment and install dependencies:

```
conda create --yes -n 20220405_chris_cerebral_malaria
mamba install -n 20220405_chris_cerebral_malaria --yes --file requirements.txt
```

## Run pipeline

Unless you have already done so, activate the environment:

```
conda activate 20220405_chris_cerebral_malaria
```

Run pipeline:

```
snakemake -p --dryrun --jobs 11 \
    --config fastqdir=/export/projects/III-data/wcmp_bioinformatics/db291g/data \
             sample_sheet=$PWD/sample_sheet.tsv \
             species=$PWD/species.tsv \
    --directory /export/projects/III-data/wcmp_bioinformatics/db291g/projects/20220405_chris_cerebral_malaria
```

Configuration option:

* `fastqdir`: Path to fastq files listed in the sample sheet.

* `sample_sheet`: Tab separated file of sample information. Edit columns
  fastq_r1/2 with the appropriate path

* `species`: Tab separated files of species to analyse (most likely, no need to edit this)

* `--directory`: Output directory

