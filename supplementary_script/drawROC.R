library(ggplot2)

## Read the data

curatedDatabase <- read.table(file = "curatedDB_dataframe", header = TRUE)
mirCureResult <- read.delim(file = "Result.txt")
curatedDatabaseID <- as.vector(curatedDatabase$ID)
mirCureResultID <- as.vector(mirCureResult$ID)

## Set a threshold and we will train MirCure based on this threshold

threshold <- 5
FPResult <- c()
TPResult <- c()
for (k in 1: 50) {
  cat(k)
  cat('\n')
  ## set to 0 after each loop
  truePosistive <- 0
  falsePositive <- 0
  
  for (i in 1:nrow(mirCureResult)) {
    find <- FALSE
    for (j in 1:nrow(curatedDatabase)) {
      
      ## if we could find the reasult in curated database, we would increase the truePositve reuslt, and swith the find to TRUE
      if (curatedDatabaseID[j] == mirCureResultID[i] & mirCureResult$Score[i] >= threshold) {
        truePosistive = truePosistive + 1
        find <- TRUE
      }
    }
    ## if we could not find the that candidate in mirGeneDB and the result is higher than the threshold, we would increase the false negative result
    if (find == FALSE & mirCureResult$Score[i] >= threshold) {
      falsePositive = falsePositive + 1
    }
  }
  threshold = threshold - 0.25
  falsePositiveRate <- falsePositive / (nrow(mirCureResult) - nrow(curatedDatabase))
  truePosistiveRate <- truePosistive / nrow (curatedDatabase)
  
  
  ## save the result
  FPResult[k] <- falsePositiveRate
  TPResult[k] <- truePosistiveRate
}

## create the data frame
dataResult <- data.frame(TPResult,FPResult)
colnames(dataResult) <- c("True Positive Rate", "False Positive Rate")


## draw the plot
ggplot(data = dataResult) +
  geom_line(aes(y = `True Positive Rate`, x= `False Positive Rate`))
