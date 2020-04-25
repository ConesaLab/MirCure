# MirCure; Quality control and curation of microRNAs

- Guillem Ylla
- Tianyuan Liu

## Installation

 - MirCure requires 2 external dependencies.
   - [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/#download).
   - And [LATEX](https://www.latex-project.org/get/) for generating df reports (optional). 



- Install MirCure as R package. 

```
library(devtools)
install_github("guillemylla/MirCureApp")
```




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
	
	
## Auxiliary scripts
 
 - MiRbase_2_MirCure.R
 - mirDeep2_2_MirCure
