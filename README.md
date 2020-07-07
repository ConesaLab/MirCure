# MirCure: A tool for quality control, filter, and curation of microRNAs of animals and plants

![GitHub last commit](https://img.shields.io/github/last-commit/ConesaLab/MirCure)
![pubstatus](https://img.shields.io/badge/publication:-Accepted-orange)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

- Guillem Ylla
- Tianyuan Liu
- Ana Conesa



### About MirCure

MirCure is a computational application to assist on filtering and curating microRNA annotations obtained from databases or  *de novo* miRNA annotation tools.

Given a list of miRNA candidates, MirCure evaluates a number of miRNA-specific features based on gene expression, biogenesis, and conservation data, and generates a score that can be used to discard poorly supported miRNA annotations. MirCure can also curate and adjust the annotation of the 5p and 3p arms based on user-provided small RNA-seq data. We evaluated MirCure on a set of manually curated animal and plant microRNAs and demonstrated great accuracy. Moreover, we show that MirCure can be used to revisit previous *bona fide* miRNAs annotations to improve microRNA databases.

We provide MirCure as a Shinny App, which can easily be executed as an R function or Docker container.


## Video Tutorial

- We make a tutorial video to help you to use MirCure:https://youtu.be/6rBHaRwqAwA

## Introduction

- The miRNAs are short (~22nts) non-coding RNAs that bind to target messenger RNA (mRNAs) to block their translation and trigger degradation. Due to their precise and powerful mechanism of action, miRNAs are involved in the regulation of a myriad of biological processes and are involved in numerous diseases.
- Prediction of miRNA genes is a common step in genome annotation pipelines. Now, there are lots of false-positives in the prediction step.
- Currently, the few sources of reliable miRNA annotations are those that went through an extensive manual curation process. This process involves running different algorithms and collecting different pieces of information that are then put together for human assessment. This is a tedious and time-consuming effort.
- To assist researchers in the selection of bona fide miRNAs out of a list of putative miRNA annotations, we developed MirCure. MirCure is a Shiny tool for quality control and filtering of both existing miRNA annotations and de novo miRNA predictions. MirCure obtains a quality score based on the automatic evaluation of different sources of miRNA evidence that can be used to accept, reject, or revisit miR-NA calls. Moreover, MirCure generates a visual report for each miRNA candidate in which all the relevant miRNA biogenesis criteria are represented, thereby greatly facilitating the curation process. We validate MirCure on a set of extensively curated miRNA annotations both for animals and plants. To our knowledge, MirCure is the first tool to perform automatic quality control of miRNA annotations based on small RNA-seq and strict biogenesis criteria.

## Basic Strategy

- MirCure uses small RNA-seq coverage, structural information, and sequence conservation to calculate a quality score that can be use to select robustly supported miRNA annotations.

- Example of some of the features that MirCure analyzes:

<p align="center">
<img src="https://github.com/ConesaLab/MirCure/blob/master/Figures/Figure3_final.png" width="450">
	</p>
<p align="center">
Figure 1
	</p>

- A real microRNA should have:
	1. Hairpin-like  secondary structure of the precursor sequence.
	2. Two miRNA between 20-26 nts long with at least 16 complementary nts.
	3. The two-arm (mature/star) arm should have two nts overhang and both 3' extremes.
	4. RNA-seq evidence supporting the mature/star arms.
	5. ow proportion of reads within the loop region and precursor flanking regions.
  6. In plants,  we expect higher seqeunce  conservation.


<p align="center">
<img src="https://github.com/ConesaLab/MirCure/blob/master/Figures/Figure1_MirCure_Pipeline_Final.jpg" width="750">
	</p>
<p align="center">
Figure 2
	</p>


## How to Install and Run MirCure

There are two main options to run MirCure:

1.  Run MirCure Docker image (**recommended**)
2.  Install and run MirCure as  *R package*.

### Run  MirCure  Docker Image

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/conesalab/mircure)


- You just need to install the [Docker engine](https://docs.docker.com/engine/install/).
- Run MirCure with the following command:

```
docker run --rm -p 3838:3838 \
-v LOCAL_Genome_DIR:/srv/shiny-server/data/genomes \
-v LOCAL_BAMFILE_DIR:/srv/shiny-server/data/bamfiles \
-v LOCAL_DIR_TO_SAVE_MIRCUREREPORTS:/srv/shiny-server/reports\
    conesalab/mircure
```

- Where:
	- **LOCAL_Genome_DIR** is the path to directory where the **genome fasta file** is located.
	- **LOCAL_BAMFILE_DIR** is the path to directory where the **bam file** is located.
	- **LOCAL_DIR_TO_SAVE_MIRCUREREPORTS** path to directory where **MirCure** pdf reports will be saved.

- Open the URL that the terminal will return (typically [http://[::]:3838](http://[::]:3838)) in a **web browser**.


![](Figures/mircure_url.png)

### MirCure as R package

 - MirCure requires 2 external dependencies.

   - [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/#download).
   - And [LATEX](https://www.latex-project.org/get/) for generating pdf reports (optional).

- Install MirCure as R package.

```
library(devtools)
install_github("ConesaLab/MirCure")
```
- Launch MirCure as R function

```
library("MirCure")
runMirCure()
```

## Example Files

Download the example files to test MirCure.
https://uflorida-my.sharepoint.com/:u:/g/personal/tianyuan_liu_ufl_edu/EecwMk1SxX9EluiTkOojjHQBIBSbPBe12hA9ZTMTnXf0Dg?e=Uk4Bet


## Input Data for MirCure

- MirCure Inputs
	1. Genome (fasta).
	2. microRNA annotations to evaluate, in 3 gff files. One with the precursor cordinates, and two with the 2 miRNA arms (mature and star, or 3' and 5').
	3. Small RNA-seq data mapped to the genome (bam)

- Some Scripts to help you create a bam file in [**How to create bam files.RTF**](https://github.com/ConesaLab/MirCure/blob/master/How%20to%20create%20bam%20files.rtf)

## Auxiliary scripts

 - MiRbase_2_MirCure.R
 - mirDeep2_2_MirCure

## Citation

Ylla, G., Tianyuan, L., & Ana, C. (2020). MirCure: A tool for quality control, filter, and curation of microRNAs of animals and plants.
