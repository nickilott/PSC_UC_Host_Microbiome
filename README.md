
# Tissue-dependent transcriptional and bacterial associations in primary sclerosing cholangitis-associated inflammatory bowel disease

## Overview

This repository contains code that was used to analyse data from [pu](publication_url). Briefly, the aim of the study was to investigate changes in gene expression in PSC/UC relative to healthy controls and patients with UC. In addition, we also analysed 16S rRNA amplicon sequencing data in order to look for genus abundance changes in these same disease cohorts. I have attempted to make the analyses performed in the manuscript as reproducible as I can and in theory all of the figures from the manuscript (as well as additional analyses) should be able to be reproduced by following the instructions below.

## The scripts

The scripts that will be run to recreate the analysis are a series of Rmarkdown files that you can find in the Rmd/ directory of this repository. The PSC_UC_Host_Microbiome.Rmd is the main file that pulls together the others from Rmd/ into a single output report.


## Installation

It should be straightforward to re-run the analyses using the same environment as was used in the primary analyses as I have set up a conda environment .yml that will recreate my environment with all of the neccessary dependencies. Therefore to re-create this make sure you have conda installed, clone this repository and do:


```
cd <path-to-repo>/PSC_UC_Host_Microbiome
conda env create -f envs/psc_uc_host_microbiome_enviroment.yml
conda activate psc_uc_host_microbiome
```

This will create an environment with all of the relevant dependencies installed.

## Install CombatSeq

I have to admit that I had issues when trying to install sva-devel (that contains CombatSeq) from github. Therefore I did a little workaround that was suggested on github issues. This suggestion was to download the code and run the individual R scripts without installing. This worked for me and so in order to re-run my analyses you will have to do the same. From inside the PSC_UC_Host_Microbiome/ directory, clone [sva-devel](https://github.com/jtleek/sva-devel). Having done this the scripts that are required by my scripts will be in the correct location.

```
git clone https://github.com/jtleek/sva-devel.git
```

## Download metadata from figshare

In order to run the analysis you will also have to have access to patient metadata, counts tables and annotation files. Download the metadata.zip, counts.zip and annotations.zip from the project at [figshare](https://figshare.com/projects/Tissue-dependent_transcriptional_and_bacterial_associations_in_primary_sclerosing_cholangitis-associated_inflammatory_bowel_disease/97364). These should be downloaded into the PSC_UC_Host_Microbiome/ directory and unzipped there. 


## Pre-made html report

While I fully expect you to be able to re-run the analysis (maybe), I have also put the final report in this repository (report/PSC_UC_Host_Microbiome.html) so that you can browse the analyses that were undertaken. This report is split into sections that largely correspond to the order of the manuscript although there are additional analyses in here that did not make it into the manuscript but are kept for reasons of transparency.


## Run the analyses

Once you have activate the correct conda environment, downloaded sva-devel and downloaded the relevant files from figshare make sure you are in the repository directory and you should simply be able to start R and run:

```
rmarkdown::render("Rmd/PSC_UC_Host_Microbiome.Rmd", output_format="html_document")
```

It may take some time to run but if successful you will find the report written out as Rmd/PSC_UC_Host_Microbiome.html.

