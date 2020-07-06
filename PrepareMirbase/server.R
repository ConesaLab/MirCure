# Define server logic to summarize and view selected dataset ----
library(DT)
library(rtracklayer)
library(IRanges)

###in cluster
#system("module load viennarna")

memory.size(max=TRUE)

options(shiny.maxRequestSize = 10000*1024^2)



server <- function(input, output, session) {

  output$visualize <- renderTable({
    shiny::validate(
      need(!is.null(input$inputfile1$datapath), message="Input miRBase gff3 (ftp://mirbase.org/pub/mirbase/CURRENT/genomes)")
    )
    mirbasegff3 <<- readGFF(input$inputfile1$datapath)
    return(head(mirbasegff3))

  })  #visualize function ends



  observeEvent(input$run, {

    ## Get all the precursors
    #mirbase_precs<-mirbasegff3[mirbasegff3$type=="miRNA_primary_transcript",]

    # get all the miRNAS (fopr some we have mature + star, for others only mature)
    # some have 5p/3p but not all
    mirbasegff3asranges<-makeGRangesFromDataFrame(mirbasegff3, keep.extra.columns=TRUE)
    mirbasegff3asranges_precs<-mirbasegff3asranges[mirbasegff3asranges$type=="miRNA_primary_transcript",]
    mirbasegff3asranges_No_precs<-mirbasegff3asranges[!mirbasegff3asranges$type=="miRNA_primary_transcript",]

    if(max(countOverlaps(mirbasegff3asranges_precs, mirbasegff3asranges_precs, ignore.strand=FALSE))>1) {print("Error!! overlaping precursors!!!!!")}else{"No overlapping precursors!"}

    overlapsprecs_norpecs<-findOverlaps(mirbasegff3asranges_precs, mirbasegff3asranges_No_precs)
    overlapsprecs_norpecs


    precs_vector=NULL
    fiveP_vector=NULL
    threeP_vector=NULL


    withProgress(message = 'Progress', value = 0, {

      n=length(mirbasegff3asranges_precs)

      for(i in 1:n){# for each precursor

        # print(i)
        incProgress(1/n)
        prec_products<- mirbasegff3asranges_No_precs[ subjectHits(findOverlaps(mirbasegff3asranges_precs[i], mirbasegff3asranges_No_precs)), ]

        if (length(prec_products) < 3) {
        preccoordinates<-mirbasegff3asranges_precs[i]
        precs_vector<-rbind(precs_vector,as.data.frame(preccoordinates))
        }

        if(length(prec_products)==2){# if 2 arms annotated
          if(as.character(strand(prec_products[1]))=="+") {#if + strand
            fiveP_vector=rbind(fiveP_vector, as.data.frame(prec_products[ start(prec_products)==min(start(prec_products))])) ## 5p is that one with lower start value
            threeP_vector=rbind(threeP_vector, as.data.frame(prec_products[ !start(prec_products)==min(start(prec_products))]))## 3p is the other
          }else if (as.character(strand(prec_products[1]))=="-") {#if - strand
            fiveP_vector=rbind(fiveP_vector, as.data.frame(prec_products[ end(prec_products)==max(end(prec_products))])) ## 5p is that one with lower start value
            threeP_vector=rbind(threeP_vector, as.data.frame(prec_products[ !end(prec_products)==max(end(prec_products))]))## 3p is the other
          }

        } else if (length(prec_products)==1) {  ## if Only 1 arm, I need to dientify if its 5p or 3p and rpedict the other
          if(as.character(strand(prec_products[1]))=="+") {#if + strand
            if(start(prec_products)-start(mirbasegff3asranges_precs[i])  <  end(mirbasegff3asranges_precs[i])-end(prec_products)   ) {# if product closer to start=5p
              print(paste("case1",i))

              prec_products1=prec_products
              prec_products1$Name=paste0(prec_products$Name,"_5p")#add 5p to name
              fiveP_vector=rbind(fiveP_vector, as.data.frame(prec_products1)) ## 5p is that one with lower start value
              ##predict 3p:
              flanking=start(prec_products)-start(preccoordinates)
              productlength=width(prec_products)
              threeP=preccoordinates
              start(threeP)=end(preccoordinates)-flanking-productlength+1
              end(threeP)=end(preccoordinates)-flanking
              threeP$source="predicted"
              threeP$Alias="NA"
              threeP$type="miRNA"
              threeP$ID="NA"
              threeP$Name=paste0(preccoordinates$Name,"_3p")
              threeP_vector=rbind(threeP_vector, as.data.frame(threeP)) ##

            }else{
              print(paste("case2",i))
              prec_products1=prec_products
              prec_products1$Name=paste0(prec_products1$Name,"_3p")#add 3p  to name
              threeP_vector=rbind(threeP_vector, as.data.frame(prec_products1)) ## 5p is that one with lower start value
              #predict 5p:

              flanking=end(preccoordinates)-end(prec_products)
              productlength=width(prec_products)
              fiveP=preccoordinates
              start(fiveP)=start(preccoordinates)+flanking
              end(fiveP)=start(preccoordinates)+flanking+productlength-1
              fiveP$source="predicted"
              fiveP$Alias="NA"
              fiveP$type="miRNA"
              fiveP$ID="NA"
              fiveP$Name=paste0(preccoordinates$Name,"_5p")
              fiveP_vector=rbind(fiveP_vector, as.data.frame(fiveP)) ##
            }
          }else if (as.character(strand(prec_products[1]))=="-") {
            if(start(prec_products)-start(mirbasegff3asranges_precs[i])  >  end(mirbasegff3asranges_precs[i])-end(prec_products)   ) {# if product closer to start=5p
              print(paste("case3",i))
              prec_products1=prec_products
              prec_products1$Name=paste0(prec_products1$Name,"_5p")#add 5p  to name
              fiveP_vector=rbind(fiveP_vector, as.data.frame(prec_products1)) ## 5p is that one with lower start value
              ##predict 3p:
              flanking=end(preccoordinates)-end(prec_products)
              productlength=width(prec_products)
              threeP=preccoordinates
              start(threeP)=start(preccoordinates)+flanking
              end(threeP)=start(preccoordinates)+flanking+productlength-1
              threeP$source="predicted"
              threeP$Alias="NA"
              threeP$type="miRNA"
              threeP$ID="NA"
              threeP$Name=paste0(preccoordinates$Name,"_5p")
              threeP_vector=rbind(threeP_vector, as.data.frame(threeP)) ##
            }else{
              print(paste("case4",i))
              prec_products1=prec_products
              prec_products1$Name=paste0(prec_products1$Name,"_3p")#add 5p  to name
              threeP_vector=rbind(threeP_vector, as.data.frame(prec_products1)) ## 5p is that one with lower start value
              #predict 5p:
              flanking=start(prec_products)-start(preccoordinates)
              productlength=width(prec_products)
              fiveP=preccoordinates
              start(fiveP)=end(preccoordinates)-flanking-productlength+1
              end(fiveP)=end(preccoordinates)-flanking
              fiveP$source="predicted"
              fiveP$Alias="NA"
              fiveP$type="miRNA"
              fiveP$ID="NA"
              fiveP$Name=paste0(preccoordinates$Name,"_3p")
              fiveP_vector=rbind(fiveP_vector, as.data.frame(fiveP)) ##

            }
          }
        } else {
          next
        }
      }
    })#clsoe indocator

    ###### lets avoid repeated names adding _2,_3 etc

    Newnames= make.names(precs_vector$Name,unique=T)
    Newnames=gsub("\\.", "-", Newnames)
    precs_vector$Name<-Newnames
    fiveP_vector$Name<-paste0(Newnames,"_5p")
    threeP_vector$Name<-paste0(Newnames,"_3p")


    output$precs <- renderTable({head(precs_vector)})
    output$fiveP <- renderTable({head(fiveP_vector)})
    output$threeP <- renderTable({head(threeP_vector)})

    ## make them global to be downloeaded
    precs_vector<<-cbind(precs_vector[,1],".","miRNA_primary_transcript", precs_vector[,2:3], ".",precs_vector$strand,".", paste("ID=",precs_vector[,12],";Alias=",precs_vector[,11],";OldId=",precs_vector[,10] ,sep=""))
    fiveP_vector<<-cbind(fiveP_vector[,1],".","miRNA", fiveP_vector[,2:3], ".",fiveP_vector$strand,".", paste("ID=",fiveP_vector[,12],";Alias=",fiveP_vector[,11],";OldId=",fiveP_vector[,10] ,sep=""))
    threeP_vector<<-cbind(threeP_vector[,1],".","miRNA", threeP_vector[,2:3], ".",threeP_vector$strand,".", paste("ID=",threeP_vector[,12],";Alias=",threeP_vector[,11],";OldId=",threeP_vector[,10] ,sep=""))
})





  ###############
  ############### Download data
  # Downloadable csv of SELECTED ROWS ----
  output$downloadprecs <- downloadHandler(
    filename = "Precursors.gff3",
    content = function(file) {
      write.table(precs_vector, file, row.names = FALSE, col.names = FALSE,quote = FALSE, sep="\t")
    }
  )
  output$downloadp5p <- downloadHandler(
    filename = "5PmiRNA.gff3",
    content = function(file) {
      write.table(fiveP_vector, file, row.names = FALSE, col.names = FALSE,quote = FALSE, sep="\t")
    }
  )

  output$download3p <- downloadHandler(
    filename = "3PmiRNA.gff3",
    content = function(file) {
      write.table(threeP_vector, file, row.names = FALSE, col.names = FALSE,quote = FALSE, sep="\t")
    }
  )


}#close server
