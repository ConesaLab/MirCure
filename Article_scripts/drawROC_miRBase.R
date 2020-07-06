library(ggplot2)

curatedDatabase <- read.table(file = "curatedDB_dataframe", header = TRUE)
mirCureResult <- read.delim(file = "Resulth.txt")
curatedDatabaseID <- as.vector(curatedDatabase$ID)
mirCureResultID <- as.vector(mirCureResult$ID)



threshold <- 5
FPResult <- c()
TPResult <- c()
expression_Threshold <- 2.5
strcuture_Threshold <- 2
threshold_vector <- c()
expression_Threshold_vector <- c()
structure_Threshold_vector <- c()
for (k in 1: 80) {
  cat(k)
  cat('\n')
  ## set to 0 after each loop
  truePosistive <- 0
  falsePositive <- 0
  
  for (i in 1:nrow(mirCureResult)) {
    find <- FALSE
    for (j in 1:nrow(curatedDatabase)) {
      
      ## if we could find the reasult in curated database, we would increase the truePositve reuslt, and swith the find to TRUE
      if (curatedDatabaseID[j] == mirCureResultID[i] & (mirCureResult$Score[i] >= threshold | (mirCureResult$Score_expression_animal[i] >= expression_Threshold & mirCureResult$overhangs_score_animal[i] >= strcuture_Threshold))) {
        truePosistive = truePosistive + 1
        find <- TRUE
      }
    }
    ## if we could not find the that candidate in mirGeneDB and the result is higher than the threshold, we would increase the false negative result
    if (find == FALSE & (mirCureResult$Score[i] >= threshold | (mirCureResult$Score_expression_animal[i] >= expression_Threshold & mirCureResult$overhangs_score_animal[i] >= strcuture_Threshold) )) {
      falsePositive = falsePositive + 1
    }
  }
  expression_Threshold_vector [k] <- expression_Threshold
  structure_Threshold_vector[k] <- strcuture_Threshold
  threshold_vector[k] <- threshold
  
  threshold = threshold - 0.25
  expression_Threshold <- expression_Threshold - 0.15
  strcuture_Threshold <- strcuture_Threshold - 0.15
  falsePositiveRate <- falsePositive / (nrow(mirCureResult) - nrow(curatedDatabase))
  truePosistiveRate <- truePosistive / nrow (curatedDatabase)
  
  
  ## save the result
  FPResult[k] <- falsePositiveRate
  TPResult[k] <- truePosistiveRate
}



dataResult <- data.frame(TPResult, FPResult)
colnames(dataResult) <- c ("TPR", "FPR")
## create the data frame
ggplot(data = dataResult) +
  geom_line(aes(y = `True Positive Rate`, x= `False Positive Rate`), size = 2) +
  labs(title = "Human miRBase filter ROC curve") +
  theme(
    axis.text.x = element_text(size = 25, family = "sans"),
    axis.text.y = element_text(size = 25,family = "sans"),
    axis.title = element_text(size = 25, family = "sans"),
    plot.title = element_text(color = "black", size = 20, face = "bold",hjust = 0.5, family = "sans"),
    plot.caption = element_text(color = "black", size = 20, family = "sans"),
    legend.position = 'none'
  )



dataResultAll <-  data.frame(TPResult, FPResult, threshold_vector, structure_Threshold_vector, expression_Threshold_vector)


write.table(dataResultAll, file =  "human_miRBase_Result.txt")
