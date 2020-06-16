csv <- read.csv("mirdeep_result_mmu.csv")# read the mirDeep result
bed <- read.delim2("mouse.bed", header = FALSE)n# read the mirDeep result

# create some vector to save the information for the gff file
mature_end <- c()
mature_start <- c()
star_end <- c()
star_start<- c()
precursor_start <- c()
precursor_end <- c()
# use loop to find the position of the star sequence and the mature sequence
for (i in 1:nrow(csv)) {
  # grep the mature and star sequence in the precursor.
  mature <- gregexpr(csv$consensus.mature.sequence[i],  csv$consensus.precursor.sequence[i])
  star <- gregexpr(csv$consensus.star.sequence[i], csv$consensus.precursor.sequence[i])
  # find the sequence long of the mature and star sequence.
  mature_long <- nchar(as.character(csv$consensus.mature.sequence[i]))
  mature_long
  star_long <- nchar(as.character(csv$consensus.star.sequence[i]))
  star_long

  mature_1 <- mature[[1]][1]
  star_1 <- star[[1]][1]

  # find the start and the end position of the sequence in the genome;
  # Because we use GFF file, we need to plus one for the start;
  precursor_start[i] <- bed$V2[i] + 1
  precursor_end[i] <- bed$V3[i]
  
  start_mature <- bed$V2[i] + mature_1 - 1 + 1
  mature_start[i] <- start_mature
  end_mature <- start_mature + mature_long - 1
  mature_end[i] <- end_mature

  start_star <- bed$V2[i] + star_1 - 1 + 1
  star_start[i] <- start_star
  end_star <- start_star + star_long - 1 
  star_end[i] <- end_star
}


# create some important characters in the GFF3 file.
sequence_ID <- as.character(bed$V1)
strand <- as.character(bed$V6)
type1 <- c()
type2 <- c()
score <- c()
source <- c()
score <- c()
phase <- c()
for (i in 1:nrow(csv)){
  type1[i] <- c("miRNA")
  type2[i] <- c("miRNA_primary_transcript")
  score[i] <- c(".")
  source <- c(".")
  phase <- c(".")
}

# there are some novel miRNA, so I name them by Novel-miRNA_i
# other mature miRNA are named by their ID
ID <- c()
for (i in 1:nrow(csv)){
  if (csv$miRBase.miRNA[i] == "-") {
    ID[i] <- paste("Novel-miRNA_", i, sep = "")
  } else {
    ID[i] <- paste (as.character(csv$miRBase.miRNA[i]), sep = "")
  }
}


# create the mature GFF3 file
csv$miRBase.miRNA <- ID
matureGFF3 <- cbind(sequence_ID, source, type1, mature_start, 
                    mature_end, score, strand, phase, ID)

write.table(matureGFF3, file = "matureGFF3.gff3", sep = "\t",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

# create the star GFF3 file
starGFF3 <- cbind(sequence_ID, source, type1, star_start,
                 star_end, score, strand, phase, ID)

write.table(starGFF3, file = "starGFF3.gff3", sep = "\t",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

#create the precursor GFF3 file
precursorGFF3 <- cbind(sequence_ID, source, type2, precursor_start,
                 precursor_end, score, strand, phase, ID)


write.table(precursorGFF3, file = "precursorGFF3.gff3", sep = "\t",
            row.names = FALSE, quote = FALSE, col.names = FALSE)


mirQCApp::runMirPlot()

