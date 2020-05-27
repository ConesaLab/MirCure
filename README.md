# MirCure: A tool for quality control, filter, and curation of microRNAs of animals and plants 

- Guillem Ylla
- Tianyuan Liu
- Ana Conesa

## tutorial Video

- We make a tutoral video to help you to use MirCure:https://youtu.be/j5FmuxpZjnY

## Introduction
- The miRNAs are short (~22nts) non-coding RNAs that bind to target messenger RNA (mRNAs) to block their translation and trigger degradation. Due to their precise and powerful mechanism of action, miRNAs are involved in the regulation of a myriad of biological processes and are involved in numerous diseases. 
- Prediction of miRNA genes is a common step in genome annotation pipelines. Now, there are lots of false-positives in the prediction step.
Currently, the few sources of reliable miRNA annotations are those that went through an extensive manual curation process. This process involves running different algorithms and collecting different pieces of information that are then put together for human assessment. This is a tedious and time-consuming effort. 
- To assist researchers in the selection of bona fide miRNAs out of a list of putative miRNA annotations, we developed Mercure. Mercure is a Shiny tool for quality control and filtering of both existing miRNA annotations and de novo miRNA predictions. Mercure obtains a quality score based on the automatic evaluation of different sources of miRNA evidence that can be used to accept, reject, or revisit miR-NA calls. Moreover, MirCure generates a visual report for each miRNA candidate in which all the relevant miRNA biogenesis criteria are represented, thereby greatly facilitating the curation process. We validate MirCure on a set of extensively curated miRNA annotations both for animals and plants. To our knowledge, Mercure is the first tool to perform automatic quality control of miRNA annotations based on small RNA-seq and strict biogenesis criteria.

## Basic Strategy
<p align="center">
<img src="https://github.com/ConesaLab/MirCure/blob/master/Figures/Figure3_final.png" width="450">
	</p>
<p align="center">
Figure 1
	</p>

- A real micro-RNA should have a good structure/expression evidence shown in Figure 1:
	1. The length of both arms should between 20-26 nts and have at least 16 complementary nts.
	2. The two-arm (mature/star) arm should have two nts overhang.
	3. Enough read evidence to support mature/star arm and fewer reads in the loop/flanking region.
	
- A real micro-RNA should highly conserved in plant
- Based on our strategy, we give different scores for different miRNA candidates. Only the candidate passed the threshold could be considered as real miRNA.

## Method
- The MirCure toolset basically consists of three software components, the MirCure R application, and two supporting scripts to convert source files to the MirCure required formats. MirCure is written in R with a Shiny interface and runs both in desktop com-puters and HPC clusters. To facilitate the installation process, Mir-Cure is also available as an R package that automatically installs all R dependencies . MirCure requires miRNA genome annotation data provided as three gff3 files, describing precursor annotations and annotations for the two miRNA arms. The gff3 files corresponding to the miRNA arms may use 5’/3’ or mature/star as an-notation styles. When arms are annotated as 3’/5’, MirCure con-siders as mature the arm with the highest number of small RNA-seq reads. Additionally, the PreapareMirbase.R script is used to predict miRNA arm annotations from mirRBase data, where this information is frequently missing. Finally, we include mirDeep_2_mirCure.R to format the output of the mirDeep2 to the MirCure format. MirCure also requires as input files the genome sequence (fasta) of the organism and a bam file with small RNA-seq data mapped to this genome sequence. 
- The MirCure software runs the pipeline described in Figure 2, that implements state-of-the-art guidelines for robust miRNAs annotation in animals and plants, by combining expression evidence (i.e. small RNA-seq data) and biogenesis information. This pipeline first adjusts miRNA annotations based on small RNA-seq data and then evaluates three different aspects of miRNA quality: secondary structure, gene expression, and conservation. For each of these aspects, a score is calculated, and relevant graphical outputs are generated. The final quality scores assigned to each miRNA candidate is a function of different sub-scores calculated in the three MirCure steps . The weight of each subscore on the final score calculation is different for animals and plants due to the differences in their biogenesis.

<p align="center">
<img src="https://github.com/ConesaLab/MirCure/blob/master/Figures/Figure1_MirCure_Pipeline_Final.jpg" width="750">
	</p>
<p align="center">
Figure 2
	</p>


## Installation

 - MirCure requires 2 external dependencies.
   - [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/#download).
   - And [LATEX](https://www.latex-project.org/get/) for generating df reports (optional). 

- Install MirCure as R package. 

```
library(devtools)
install_github("ConesaLab/MirCure")
```
## Example files

Download the example files to run and understand how to use MirCure.
https://uflorida-my.sharepoint.com/:u:/g/personal/tianyuan_liu_ufl_edu/EV9JbkyGvqxPj1-RZHFSltUBG8M4Hg0K9UMMBnsN32U7UA?e=ZGvTVg

## Run MirCure

- Launch MirCure user-graphic interface 

```
library("MirCure")
runMirCure()
```

- MirCure Inputs
	1. Genome (fasta).
	2. microRNA annotations to evaluate, in 3 gff files. One with the precursor cordinates, and two with the 2 miRNA arms (mature and star, or 3' and 5'). 
	3. Small RNA-seq data mapped to the genome (bam)

- Some Scripts to help you create a bam file in **How to create bam files.RTF**
## Auxiliary scripts
 
 - MiRbase_2_MirCure.R
 - mirDeep2_2_MirCure
