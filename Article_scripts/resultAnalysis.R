setwd("~/Desktop/human")
library(rPraat)
library(testit)
library(rtracklayer)
library(GenomicAlignments)

pre <- readGFF("miRNApre.gff3") ##  Load the precursor file 
csv <- read.csv(file = "FilteredMiRNAs.csv" ) # Load the filtered miRNA file
curatedDB <- readGFF(file = "hsa.gff") # read the data from the curated database(MirGeneDB)

chr <- curatedDB$seqid
chr <- gsub("chr", "",chr) # remove "chr" since there is no chr in the miRNApre.gff3
curatedDB$seqid <- chr

curatedDBpre <- c()
curatedDB5P <- c()
curatedDB3P <- c()

# extract the 5p/3p/pre annotation from the curated database
for (i in 1 : nrow(curatedDB)) {
  curatedDBname <- strsplit (curatedDB$ID, "_")[[i]][2]
  curatedDBname <- substr(curatedDBname, start = 1, stop = 1)
  if (curatedDBname == "p") {
    curatedDBpre[i] = i
  } else if (curatedDBname == "5") {
    curatedDB5P[i] = i
  } else {
    curatedDB3P[i] = i
  }
}

curatedDBpre <- curatedDBpre[!is.na(curatedDBpre)]
curatedDB5P <- curatedDB5P[!is.na(curatedDB5P)]
curatedDB3P <- curatedDB3P[!is.na(curatedDB3P)]


curatedDBpre <- curatedDB[curatedDBpre,]
curatedDB5P <- curatedDB[curatedDB5P,]
curatedDB3P <- curatedDB[curatedDB3P,]


################################## compare precursor ###########################################

# Use GRanges to compare the position of the mature star sequence.
miRBasePre <-as(pre, "GRanges")
curatedDBpre <- as(curatedDBpre, "GRanges")
head(miRBasePre)
#GoodmiRBasePre <- as.data.frame(intersect(miRBasePre,curatedDBpre))
#GoodmiRBasePre <- GoodmiRBasePre[GoodmiRBasePre$width >= 50,]

GoodmiRBasePre <- as.data.frame(findOverlaps(miRBasePre, curatedDBpre, type = "any", minoverlap = 50L)) # If there are at least 50 common nts between the two sequence, we think they are the same miRNA
curatedDBpre <- miRBasePre[unique(GoodmiRBasePre$queryHits),]
curatedDBpre <- as(curatedDBpre,"GRanges")

### Start to calculate true positives/ false positives/ fales negatives/ true negatives
library(rPraat)
library(testit)
library(rtracklayer)
library(GenomicAlignments)
library(eulerr)

filterResult <- read.csv("FilteredMiRNAs.csv")

ID <- as.character(filterResult$ID)
strand <- c()
chromosome <- c()
start <- c()
end <- c()

for (i in 1 : nrow(filterResult)){
  chromosome[i] <- strsplit(as.character(filterResult$Loci[i]), ":|-")[[1]][1]
  start[i] <- strsplit(as.character(filterResult$Loci[i]), ":|-")[[1]][2]
  end[i] <- strsplit(as.character(filterResult$Loci[i]), ":|-")[[1]][3]
  if (strsplit(as.character(filterResult$Loci[i]), ":|-")[[1]][4] == ""){
    strand[i] <- "-"
  } else {
    strand[i] <- strsplit(as.character(filterResult$Loci[i]), ":|-")[[1]][4]
  }
}

type <- c()
score <- c()
source <- c()
phase <- c()
for (i in 1:nrow(filterResult)){
  type[i] <- c("miRNA")
  score[i] <- c(".")
  source <- c(".")
  phase <- c(".")
}

filterGFF <- as.data.frame(cbind(chromosome, source, type, start, end, score, strand, phase, ID))
filterRange <- as(filterGFF, "GRanges")

goodFilter <- as.data.frame(findOverlaps(filterRange, curatedDBpre, type = "any", minoverlap = 50L)) 
truePositive <- filterGFF[goodFilter$queryHits,] 
falsePositive <- filterGFF[-c(goodFilter$queryHits),];falsePositive
falseNegative <- as.data.frame(curatedDBpre[-c(goodFilter$subjectHits),]);falseNegative


################################################analysis##################################################
library(VennDiagram)

draw.triple.venn(area1 =10, area2 = 20,
                 n12 = 5, category = c("miRBase", "mirQC"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"),euler.d = TRUE, scaled=TRUE,
                 cat.default.pos = "text")

######################statistical analysis############

falsePositiveRate <- round(nrow(falsePositive) / nrow(filterGFF),digits = 2)
falseNegativeRate <- round(nrow(as.data.frame(falseNegative)) / nrow(as.data.frame(curatedDBpre)), digits = 2)
# Load ggplot2
library(ggplot2)

# Create data
data <- data.frame(
  name=c("FNR","FPR") ,  
  value=c(falseNegativeRate, falsePositiveRate)
)
ggplot(data, aes(x=name, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label=data$value), vjust=0.5, hjust=1.1, color="white", size=4.5)+
  theme_minimal()







