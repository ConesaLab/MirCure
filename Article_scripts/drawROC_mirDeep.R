library(ggplot2)
library(GenomicAlignments)
miRdeep <- read.delim(file = "precursorGFF3.gff3", header = FALSE)
colnames(miRdeep) = c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "ID")
resultDelim <- read.delim(file = "Result_deep.txt")
curatedDataBase <- read.delim(file = "curatedDB_dataframe_hsa")
curatedDataBase <- as(curatedDataBase, "GRanges")
curatedDataBase1 <- as.data.frame(curatedDataBase)

threshold <- 5
expression_Threshold <- 2.5
strcuture_Threshold <- 2
FPResult <- c()
TPResult <- c()
threshold_vector <-c()
expression_Threshold_vector <- c()
strcuture_Threshold_vector <- c()

for (i in 1:80) {
  
  filterResult <- resultDelim$Score >= threshold | resultDelim$Score_expression_animal >= expression_Threshold & resultDelim$overhangs_score_animal >= strcuture_Threshold
  filterDeepGFF <- miRdeep[filterResult,]
  filterDeep <- as(filterDeepGFF, "GRanges")
  
  validation <- as.data.frame(findOverlaps(filterDeep, curatedDataBase, type= "any", minoverlap = 50L))
  
  truePosistive <- nrow (filterDeepGFF[validation$queryHits,])
  falsePositive <- nrow(filterDeepGFF[-c(validation$queryHits),])
  
  falsePositiveRate <- falsePositive / (nrow(miRdeep) - nrow(curatedDataBase1))
  truePosistiveRate <- truePosistive / nrow (curatedDataBase1)
  
  threshold_vector[i] <- threshold
  expression_Threshold_vector [i] <- expression_Threshold
  strcuture_Threshold_vector [i] <- strcuture_Threshold
  ## save the result
  FPResult[i] <- falsePositiveRate
  TPResult[i] <- truePosistiveRate
  expression_Threshold = expression_Threshold - 0.15
  strcuture_Threshold = strcuture_Threshold -0.1

  threshold = threshold - 0.25
}
dataResult <- data.frame(TPResult,FPResult)
colnames(dataResult) <- c("TPR", "FPR")

dataResultmmu <- data.frame(TPResult,FPResult)
colnames(dataResultmmu) <- c("True Positive Rate", "False Positive Rate")
dataResulhsa <- data.frame(TPResult,FPResult)
colnames(dataResulhsa) <- c("True Positive Rate", "False Positive Rate")

dataResultmmu$Group = c("Mouse")
dataResulhsa$Group = c("Human")

dataAll <-rbind(dataResulhsa, dataResultmmu)
## draw the plot
ggplot(data = dataAll) +
  geom_line(aes(y = `True Positive Rate`, x= `False Positive Rate`, color = Group), size = 2) +
  scale_color_manual(values = c("#00BFC4", "#C77CFF"))+
  labs(subtitle = "Figure4 \n B") +
  theme(
    legend.text = element_text(size = 15, family = "sans"),
    legend.title = element_text(color = "black", size = 18, face = "bold",family = "sans"),
    plot.subtitle = element_text(color = "black", size = 20, face = "bold",family = "sans"),
    legend.position = c(0.92,0.055),
    axis.text.x = element_text(size = 20, family = "sans"),
    axis.text.y = element_text(size = 20,family = "sans"),
    axis.title = element_text(size = 20, family = "sans"),
    plot.title = element_text(color = "black", size = 20, face = "bold",hjust = 0.5, family = "sans"),
    plot.caption = element_text(color = "black", size = 20, family = "sans"),
  )

hsa <- data.frame(FPResult,TPResult,threshold_vector,expression_Threshold_vector,strcuture_Threshold_vector)

write.table(hsa, file = "hsa_deep_result.txt")




