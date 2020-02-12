# Define server logic to summarize and view selected dataset ----
#install.packages("LncFinder", lib="/ufrc/conesa/guillemyllabou/R_libs")
#install.packages("shinyFiles", lib="/ufrc/conesa/guillemyllabou/R_libs")

###### In cluster, activate:
#system("module load viennarna")
#.libPaths(  "/ufrc/conesa/guillemyllabou/R_libs" )
######
library("GenomicAlignments")
library("Biostrings")
library("DT")
library("rtracklayer")
library("LncFinder")
library("RRNA")
library("stringr")
library("Rsubread")
library("msa")## for multiple alignments
library("shinyFiles")
library("ggplot2")

memory.size(max=TRUE)

options(shiny.maxRequestSize = 10000*1024^2)



server <- function(input, output, session) {
  #######
  ###### variables that are session-specific
  #vars to track progress
  values <- reactiveValues()

  values$successStep1<-FALSE
  values$successStep2<-FALSE
  values$successStep3<-FALSE
  values$successStep4<-FALSE
  values$successStep5<-FALSE
  values$successStep6<-FALSE

  ## select animal or plant
  observeEvent(input$specie, {
    if(input$specie == "Animal"){
      specie<<-"Animal"
      print(specie)}

    if(input$specie == "Plant"){
      specie<<-"Plant"
      print(specie)}
  })

  #precs file
  output$precs <- renderTable({
       req(input$precs)
       precsdf <- readGFF(input$precs$datapath)
       #precsdf <- read.table("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/Zma_cons_precursor.gff3")
       #precsdf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/Zma_cons_precursor.gff3")
       #precsdf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/bger/Conserved_Precursors.gff3")
       precsdf <- readGFF("/home/guillem/Documents/mirQCApp/mousedata/precursorGFF3conflict.gff3")

       return(head(precsdf))
  })

       #mat file
  output$mature <- renderTable({
    req(input$mature)
    matdf <- readGFF(input$mature$datapath)
    #matdf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/Zma_cons_mature.gff3")
    #matdf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/bger/Bger_matures.gff3")
    matdf <- readGFF("/home/guillem/Documents/mirQCApp/mousedata/matureGFF3conflict.gff3")
    return(head(matdf))
   })

         #star file
    output$star <- renderTable({
        req(input$star)
        stardf <- readGFF(input$star$datapath)
        #stardf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/Zma_cons_star.gff3")
        #stardf <- readGFF("/home/guillemyllabou/Documents/mirPlot_Shiny/v0/data/bger/Bger_stars.gff3")
        stardf <- readGFF("/home/guillem/Documents/mirQCApp/mousedata/starGFF3conflict.gff3")
        return(head(stardf))
 })

    ## user uplaoded mature/star or 5p/3p

    observeEvent(input$matureorarm, {
      if(input$matureorarm == "maturestar"){
        maturestar=TRUE
        print("maturestar")}

      if(input$matureorarm == "arm5p3p"){
        maturestar=FALSE
        print("arm5p3p")}
    })

##############################################
################# Get sequences      ##########
##############################################

### select genome

  ##if uploaded
observeEvent(input$genome0, {
  values$genomefile<<-input$genome0$datapath
  print(values$genomefile)
})
##if selected
observeEvent(input$genome1, {
    values$genomefile<<-paste("data/genomes",input$genome1, sep="/")
    print(values$genomefile)
  })
###


observeEvent(input$ButtonSeqs, {

    ## must load again the data
  ################ Checks if all files are available, avoids crash!
  shiny::validate(
    #need(input$genome0!=""  , message="missing genome"),
    need(values$genomefile!=""  , message="missing genome"),
    errorClass =  showNotification("Missing genome file", type= "error")
  )
  shiny::validate(
    need(!is.null(input$precs$datapath), message = ('Missing Precursor annotation file' )),
    errorClass =  showNotification("Missing Precursor annotation file", type= "error")
  )
  shiny::validate(
    need(!is.null(input$mature$datapath), message = ('Missing Mature annotation file')),
    errorClass =  showNotification("Missing Mature annotation file", type= "error")
  )
  shiny::validate(
    need(!is.null(input$star$datapath), message = ('Missing Star annotation file')),
    errorClass =  showNotification("Missing Star annotation file", type= "error")
  )


      ###################
  ###Processing....

  #genomefasta <- FaFile(paste("data/genomes","zma.AGPv4.full.fasta", sep="/"))
  #genomefasta <- FaFile(paste("data/genomes","Bgermanica.scaffolds.fa", sep="/"))
  genomefasta <- FaFile(paste("data/genomes","mmu.fa", sep="/"))

  genomefasta <- FaFile(values$genomefile)

  precsdf <<- try(readGFF(input$precs$datapath) )
  matdf <<- try(readGFF(input$mature$datapath) )
  stardf <<- try(readGFF(input$star$datapath) )



  #if any try fails:
  shiny::validate(
    need(class(precsdf)!="try-error", message = ('Prec gff3 error')),
    errorClass =  showNotification(" invalid gff3 file" ,type= "error")  )
  shiny::validate(
    need(class(matdf)!="try-error", message = ('Mat gff3 error')),
    errorClass =  showNotification(" invalid gff3 file" ,type= "error")  )
  shiny::validate(
    need(class(stardf)!="try-error", message = ('Star gff3 error')),
    errorClass =  showNotification(" invalid gff3 file" ,type= "error")  )



  #### Check if is is possible to get sequences form the gff3+genome files
  precseqs <<- try(getSeq(genomefasta,GRanges(precsdf$seqid,IRanges(start=as.numeric(precsdf$start ),end=as.numeric(precsdf$end)),strand =precsdf$strand ) ))
  matseqs <<- try(getSeq(genomefasta,GRanges(matdf$seqid,IRanges(start=as.numeric(matdf$start),end=as.numeric(matdf$end)),strand =precsdf$strand )))
  starseqs <<- try(getSeq(genomefasta,GRanges(stardf$seqid,IRanges(start=as.numeric(stardf$start),end=as.numeric(stardf$end)),strand =precsdf$strand )))

  #if any try fails
  shiny::validate(
    need(class(precseqs)!="try-error", message = ('Prec genome fail')),
    errorClass =  showNotification("Precursor coordinates not in genome", type= "error")  )
  shiny::validate(
    need(class(matseqs)!="try-error", message = ('Mat genome fail')),
    errorClass =  showNotification("Mature coordinates not in genome", type= "error")  )
  shiny::validate(
    need(class(starseqs)!="try-error", message = ('Star genome  fail')),
    errorClass =  showNotification("Star coordinates not in genome", type= "error")  )

 #matseqs <<- getSeq(genomefasta,GRanges(matdf$seqid,IRanges(start=as.numeric(matdf$start),end=as.numeric(matdf$end)),strand =precsdf$strand ))
  #starseqs <<- getSeq(genomefasta,GRanges(stardf$seqid,IRanges(start=as.numeric(stardf$start),end=as.numeric(stardf$end)),strand =precsdf$strand ) )


  #extrabases<<-4
  extrabases<<-input$extrabases

  precsdf_adjusted<-precsdf



  ##### Adjust precursors: all must star with the 1st nt of the 5p miRNA and end in the alst nt of the 3p miRNA
  for(i in 1:length(precseqs)){


      allcords=c( matdf$start[i], matdf$end[i], stardf$start[i], stardf$end[i])
      precsdf_adjusted$start[i]=min(allcords)
      precsdf_adjusted$end[i]=max(allcords)

  }

  precsdf_adjusted<<-precsdf_adjusted
  #precocoordsadjusted=GRanges(precsdf_adjusted$seqid,IRanges(start=as.numeric(precsdf_adjusted$start ),end=as.numeric(precsdf_adjusted$end)))
  ## we crop the precursor from mature to end, and later we add +- extrabases
  precseqs_extended <<- getSeq(genomefasta,  GRanges(precsdf_adjusted$seqid,IRanges(start=as.numeric(precsdf_adjusted$start ),end=as.numeric(precsdf_adjusted$end) ),strand =precsdf_adjusted$strand )+ extrabases )




  #output$mirnaSeqs <-renderTable({
  #mirnadf<<-data.frame(ID=as.character(precsdf$ID), mature=as.character(matseqs),star=as.character(starseqs),"Loci"=paste(precsdf$seqid,":",precsdf$start,"-",precsdf$end,":",precsdf$strand, sep=''),precursor=as.character(precseqs),precseqs_extended=as.character(precseqs_extended) )

  mirnadf<<-data.frame(ID=as.character(precsdf$ID), mature=as.character(matseqs),star=as.character(starseqs),"Loci"=paste(precsdf$seqid,":",precsdf$start,"-",precsdf$end,":",precsdf$strand, sep=''),precursor=as.character(precseqs),precseqs_extended=as.character(precseqs_extended) )

  mirnadf_toshow1<-data.frame(ID=as.character(precsdf$ID), mature=as.character(matseqs),'length mature'=width(matseqs),star=as.character(starseqs),'length mature'=width(starseqs),"Loci"=paste(precsdf$seqid,":",precsdf$start,"-",precsdf$end,":",precsdf$strand, sep=''),precursor=as.character(precseqs),precseqs_extended=as.character(precseqs_extended) )

  if (input$matureorarm == "maturestar"){
    output$mirnaSeqs <- DT::renderDataTable({ mirnadf_toshow1}  )
          }else{
            mirnadfrenamed<-mirnadf_toshow1
            colnames(mirnadfrenamed)[2]<-'5p Arm'
            colnames(mirnadfrenamed)[4]<-'3p Arm'
            output$mirnaSeqs <- DT::renderDataTable({ mirnadfrenamed}  )

              }


  values$successStep1<-TRUE
  showNotification("Done. Check Seq. Info Tab", type= "message")


  ########### penalty length
  penaltyscore1<-ifelse(width(matseqs)>23, -2, 0 )
  penaltyscore2<-ifelse(width(starseqs)>23, -2, 0 )
  penaltyscorelength_animal<<-penaltyscore2+penaltyscore1

  penaltyscore3<-ifelse(width(matseqs)>23, -2, 0 )
  penaltyscore4<-ifelse(width(starseqs)>23, -2, 0 )
  penaltyscorelength_plant<<-penaltyscore3+penaltyscore4
  ########


})## clsose button seqs

##############################################
################# Fold rpecursors   ##########
##############################################

observeEvent(input$ButtonFold, {

  ################ Checks if previous step was succesdul
  shiny::validate(
    need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
    errorClass =  showNotification("Missing succesful Step 1", type= "error")
  )


      #### requires Vienna, downalaod from https://www.tbi.univie.ac.at/RNA/index.html#download
      ########## ./configure  make , sudo make install
  withProgress(message = 'Folding precursors', value = 0, {
    overhang_animal<-NULL
    overhang2_animal<-NULL
    overhang_plant<-NULL
    overhang2_plant<-NULL
    perfectmatch_animal<-NULL
    semigood_animal<-NULL
    perfectmatch_plant<-NULL
    semigood_plant<-NULL

    ## sdefine cores animal / plant
    if (specie=="Animal"){
      perfectmatch_animal<-5
      semigood_animal<-3
    }else{
      perfectmatch_plant<-2.5
      semigood_plant<-2
    }

 n<-  nrow(mirnadf)
 for (i in 1:nrow(mirnadf)){
    #for (i in 1:5){
      #folded=run_RNAfold(as.character(mirnadf$precseqs_extended[i]), RNAfold.path = "RNAfold", detectCores(all.tests = FALSE, logical = TRUE))
      folded=run_RNAfold(as.character(mirnadf$precseqs_extended[i]), RNAfold.path = "RNAfold", parallel.cores= 4)#detectCores(all.tests = FALSE, logical = TRUE))
      coord=ct2coord(makeCt(folded[2,], folded[1,]))

      ################# Count Overhang within FOLDING ############################################
      ############################################################################################

      maturecord00<-gregexpr(RNAString(DNAString(mirnadf$mature[i])), folded[[1]][1])

       # some cases mature and star appear multiple times in prec
      if (length(maturecord00[[1]])==1) { # if matur e appears only once
        maturecord0<-maturecord00[[1]]
      }else{# if appears multiple times
        if(maturecord1[1]<starcord1[1]){ # if mature is 5p  get the first
              maturecord0<-maturecord00[[1]][1]
        }else{# if its 3p the last one
           maturecord0<-maturecord00[[1]][length(maturecord00[[1]])]
        }
      }
      maturecord1<-range(maturecord0,maturecord0+nchar(as.character(mirnadf$mature[i]))-1)



      starcord00<-gregexpr(RNAString(DNAString(mirnadf$star[i])), folded[[1]][1])
      # some cases mature and star appear multiple times in prec
      if (length(starcord00[[1]])==1) { # if star appears only once
        starcord0<-starcord00[[1]]
      }else{# if appears multiple times
        if(maturecord1[1]<starcord1[1]){ # if mature is 5p  get the last as star
            starcord0<-starcord00[[1]][length(starcord00[[1]])]
        }else{# if mature its 3p the first one
          starcord0<-starcord00[[1]][1]
        }
      }
      starcord1<-range(starcord0,starcord0+nchar(as.character(mirnadf$star[i]))-1)

      ##if 5P mature
      if(maturecord1[1]<starcord1[1]){# if 5p
        print(paste("mature is 5'",i))

        overlap=maturecord1[2]-starcord1[1]+1
        if(overlap<0) {#if mature and star don't overlap (as it should)
          colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1)) ), rep("Red",length(seq(maturecord1[1], maturecord1[2]))),rep("Black",length(seq(maturecord1[2]+1, starcord1[1]-1))),rep("Blue",length(seq(starcord1[1], starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
          foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

          foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
          foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
        }else{# If mature and star overlap (they shouldn't!) don't crash do :
          overlaFLAG=TRUE
          overlapmatstar=maturecord1[2]-starcord1[1]+1
          if (overlap>0){# if there is overlap

          colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1)) ), rep("Red",length(seq(maturecord1[1], maturecord1[2]-overlapmatstar))),rep("Orange",length(seq(maturecord1[2]-overlapmatstar+1, starcord1[1]+overlapmatstar-1))),rep("Blue",length(seq(starcord1[1]+overlapmatstar, starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
          foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

          foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
          foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
          }
          if (overlap==0){# if they dont overlap, but there is no loop

            colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1)) ), rep("Red",length(seq(maturecord1[1], maturecord1[2]-overlapmatstar))),rep("Blue",length(seq(starcord1[1]+overlapmatstar, starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
            foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

            foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
            foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
          }
        }



      }else{   ##################### Mature 3'
        print(paste("mature is 3'", i))

        overlap=starcord1[2]+1-maturecord1[1]

        if(overlap<0){ # if no overlap

          colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep("Blue",length(seq(starcord1[1], starcord1[2]))),rep("Black",length(seq(starcord1[2]+1, maturecord1[1]-1))),rep("Red",length(seq(maturecord1[1], maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
          foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

          foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
          foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)

        }else{# If mature and star overlap (they shouldn't!) don't crash do :

          if (overlap>0){
            colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep("Blue",length(seq(starcord1[1], starcord1[2]-overlapmatstar2))),rep("Orange",length(seq(starcord1[2]+overlapmatstar2+1, maturecord1[1]-overlapmatstar2-1))),rep("Red",length(seq(maturecord1[1]-overlapmatstar2, maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
            foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

            foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
           foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
          }
          if (overlap==0){# if overlap is zero but also no gap
            colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep("Blue",length(seq(starcord1[1], starcord1[2]))),rep("Red",length(seq(maturecord1[1], maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
            foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

            foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
            foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
          }

        }
      }




      #### select first miRNA  (mature or star) nucelotide and its complemenrtaey
        ## 1st I get the 1st nucleotide of the 1st miRNA (mature or star)
        firstmirnanuc5p<-min(which(foldingtable$color!="Black"))
        DroshaCut5<-foldingtable[(firstmirnanuc5p-2):(firstmirnanuc5p+1),]

        ### Check 1st cleaveage by Drosha
        foldingtable_2<-foldingtable
        foldingtable_2$dotnumeric<-0
        foldingtable_2[foldingtable_2$dots=="(",]$dotnumeric<-(-1)
        foldingtable_2[foldingtable_2$dots==")",]$dotnumeric<-1

        foldingtable_2[firstmirnanuc5p,]
        foldingtable_2

        #Findmatchingupstream(firstmirnanuc5p+1)

        ## find the matching base!
        Findmatchingupstream<-function(querynt){
          for(pos in 1:(nrow(foldingtable_2)-querynt)){
            if(sum(foldingtable_2$dotnumeric[(querynt+1) : (querynt+pos) ]) == 1){
              complement=pos+querynt
              return(complement)
              break()
            }
          }
          return(NULL)
        }

        ## I have the 1st nt of the miRNA located at the 5'

        if(foldingtable_2[firstmirnanuc5p,"dots"]!="."){#if 1st one, has a complementary
            Compl_to_firstmirnanuc5p <- Findmatchingupstream(firstmirnanuc5p) ## Get complementary to first

            if(foldingtable_2[Compl_to_firstmirnanuc5p,"color"]!="Black" & foldingtable_2[Compl_to_firstmirnanuc5p,"color"]!=foldingtable_2[firstmirnanuc5p,"color"] ){# If complementary is not same color not black
              if(foldingtable_2[Compl_to_firstmirnanuc5p+1,"color"]!="Black" & foldingtable_2[Compl_to_firstmirnanuc5p+2,"color"]!="Black" ){# if 2 upstream are also colored
                  if(sum(str_count(foldingtable_2[(firstmirnanuc5p-2):(firstmirnanuc5p+1),"dots"], "\\(" ))==4 ){ # if first, second and 2 previous all matched is perfect
                    print("Perfect Drosha")
                    overhang_animal<-rbind(overhang_animal, c("Perfect 5p Drosha cleavage", perfectmatch_animal) )
                    overhang_plant<-rbind(overhang_plant, c("Perfect 5p cleavage", perfectmatch_plant))
                  }else{#they will at least have 1 complementary (the first)
                    print("Acceptable Drosha")
                    overhang_animal<-rbind(overhang_animal, c("Acceptable 5p Drosha cleavage", perfectmatch_animal/0.9) )
                    overhang_plant<-rbind(overhang_plant, c("Acceptable 5p cleavage", perfectmatch_plant/0.9))
                  }
              }else if(foldingtable_2[Compl_to_firstmirnanuc5p+1,"color"]!="Black") {#If only 1 upstream colored, and 1st is paired
                  print("Not great Drosha 1")
                  overhang_animal<-rbind(overhang_animal, c("Weak 5p Drosha cleavage 1", perfectmatch_animal/0.5) )
                  overhang_plant<-rbind(overhang_plant, c("Weak 5p cleavage 1", perfectmatch_plant/0.5))
              }else{
                print("Bad Drosha 1")
                overhang_animal<-rbind(overhang_animal, c("Bad Drosha cleavage 1", 0) )
                overhang_plant<-rbind(overhang_plant, c("Bad Drosha cleavage 1" , 0))
              }# none upstream colored
            }else if(foldingtable_2[Compl_to_firstmirnanuc5p,"color"]=="Black" ) { #if complement of the first is a black
                print("still calculating")
                overhang_animal<-rbind(overhang_animal, c("Calculating 5p Drosha cleavage", perfectmatch_animal) )
                overhang_plant<-rbind(overhang_plant, c("Calculating 5p cleavage", perfectmatch_plant))
            }else{ #complement is same color super bad
              print("Bad Drosha 2")
              overhang_animal<-rbind(overhang_animal, c("Bad Drosha cleavage 2", 0) )
              overhang_plant<-rbind(overhang_plant, c("Bad Drosha cleavage 2" , 0))
            }
          }else if( foldingtable_2[firstmirnanuc5p+1,"dots"]!="."){# If  first doesnt have coplementary, check 2nd
              print("Could still be good")
              overhang_animal<-rbind(overhang_animal, c("Could still be good", 0) )
              overhang_plant<-rbind(overhang_plant, c("Could still be good" , 0))
            }else{ #first 2 ones no complement
              print("Bad Drosha 3")
              overhang_animal<-rbind(overhang_animal, c("Bad Drosha cleavage 3", 0) )
              overhang_plant<-rbind(overhang_plant, c("Bad Drosha cleavage 3" , 0))
            }





    ########### Check 2ns cleaveage (the loop) by Dicer
    # if no loop = bad cleavaege
        if (overlap>=-2){
          print("Overlap mature and star")
          overhang2_animal<-rbind(overhang2_animal,c("Bad 3p cleavage overhang (Dicer cutting) 1", 0))
          overhang2_plant<-rbind(overhang2_plant,c("Bad 3p cleavage overhang (Dicer cutting) 1", 0))
        }else{
          firstcolor<-foldingtable[min(which(foldingtable$color!="Black")),"color"]
          lastmirnanuc3p=max(which(foldingtable$color==firstcolor))

          DicerCut3<-foldingtable[(lastmirnanuc3p-1):(lastmirnanuc3p+2),]


          firstmirnanuc5p<-min(which(foldingtable$color!="Black" & foldingtable$color!=firstcolor) )
          DicerCut5<-foldingtable[(firstmirnanuc5p-2):(firstmirnanuc5p+1),]


          ## if the matching ones is black, +2 should be Blue or Red
          if( all(DicerCut5$dots==")" &  DicerCut3$dots == "(" )){ ## if all 8 are paired
            if(all(DicerCut5$closeprent[1:2] == rev(DicerCut3$openprent[1:2]) | DicerCut5$closeprent[3:4] == rev(DicerCut3$openprent[3:4]))){#if perfectly pared 2 to 2
              print("Perfect 3p extreme overhang 1")
              overhang2_animal<-rbind(overhang2_animal, c("Perfect 3p cleavage overhang (Dicer cutting) 1", perfectmatc2_animal) )
              overhang2_plant<-rbind(overhang2_plant, c("Perfect 3p cleavage overhang (Dicer cutting) 1", perfectmatch_plant))
            }else if (all(table(DicerCut5$closeprent %in% DicerCut3$openprent)[2]>=2)){#only 1 hangs
              print("Not Perfect 3p extreme overhang 1")
              overhang2_animal<-rbind(overhang2_animal, c("Not perfect 3p cleavage overhang (Dicer cutting) 1", perfectmatch_animal/2) )
              overhang2_plant<-rbind(overhang2_plant, c("Not perfect 3p cleavage overhang (Dicer cutting) 1", perfectmatch_plant/2))
              }else{print("Bad mature 3p overhang 1.2")
                overhang2_animal<-rbind(overhang2_animal,c("Bad 3p cleavage overhang (Dicer cutting) 1.2", 0))
                overhang2_plant<-rbind(overhang2_plant,c("Bad 3p cleavage overhang (Dicer cutting) 1.2", 0))

            }

          #}else if(!is.na(all(DicerCut5$closeprent[1:2] == rev(DicerCut3$openprent[1:2]) | DicerCut5$closeprent[3:4] == rev(DicerCut3$openprent[3:4])) )){
          }else if( sum(str_count(DicerCut3$dots, "\\("))>2 & sum(str_count(DicerCut5$dots, "\\)"))>2     ){ # if thelastone before the cut has pair and the next two also, is good!

             # if(all(DicerCut5$closeprent[1:2] == rev(DicerCut3$openprent[1:2]) | DicerCut5$closeprent[3:4] == rev(DicerCut3$openprent[3:4])) ){ # if thelastone before the cut has pair and the next two also, is good!
            if (all(table(DicerCut5$closeprent %in% DicerCut3$openprent)[2]>=2)){#only 1 hangs
               print("Good 3p not perfect overhang 2")
               overhang2_animal<-rbind(overhang2_animal, c("Not perfect 3p cleavage overhang (Dicer cutting) 2", perfectmatch_animal/2) )
               overhang2_plant<-rbind(overhang2_plant, c("Note perfect 3p cleavage overhang (Dicer cutting) 2", perfectmatch_plant/2))
              }else{
                print("Bad mature 3p overhang 3.1")
                overhang2_animal<-rbind(overhang2_animal,c("Bad 3p cleavage overhang (Dicer cutting) 3.1", 0))
                overhang2_plant<-rbind(overhang2_plant,c("Bad 3p cleavage overhang (Dicer cutting) 3.1", 0))

              }
            }else{
            print("Bad mature 3p overhang 2")
            overhang2_animal<-rbind(overhang2_animal,c("Bad 3p cleavage overhang (Dicer cutting) 3", 0))
            overhang2_plant<-rbind(overhang2_plant,c("Bad 3p cleavage overhang (Dicer cutting) 3", 0))

          }
        }
 ##############################################        ##############################################        ##############################################

      ###Lets adjust image parameters depending on length

        if(nchar(folded[1,]) > 180){

          jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
            par(mar=c(0.1,0.1,0,0.1))
            RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
              seqcols=c(2,5),labTF=FALSE,
              pointSize = 1, lineWd = 1, nt=T,
              dp=1, tsize=0.5)
            dev.off()
          }else if( nchar(folded[1,]) > 100){
            jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
              par(mar=c(0.1,0.1,0,0.1))
              RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
                seqcols=c(2,5),labTF=FALSE,
                pointSize = 2, lineWd = 1, nt=T,
                dp=1, tsize=0.8)
                dev.off()
          }else{
              jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
                par(mar=c(0.1,0.1,0,0.1))
                RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
                seqcols=c(2,5),labTF=FALSE,
                pointSize = 2.3, lineWd = 1, nt=T,
                dp=1, tsize=1)
                dev.off()
          }

   incProgress(1/n, detail = paste("Prec", i, "of", n))
    }
  })

  print(paste("Good here!! 1", i))
  foldingFigs<-paste("<img src=\"images/",mirnadf$ID,"_fold.jpg\" width=\"500\ height=\"180\"></img>", sep="")



  ## scores animal / plant
  if (specie=="Animal"){
    print("Making dataframe animal")
    print(overhang_animal)
    print(overhang2_animal)
    mirnadf_folding<-cbind(mirnadf,foldingFigs ,"5'Cleavage"=overhang_animal[,1],"3'Cleavage"=overhang2_animal[,1] )
    output$mirnaSeqswithplots <-  DT::renderDataTable({ mirnadf_folding[,c(1,2,3,7,8,9)]},  escape = FALSE )
    overhangs_score_animal<<- as.numeric(overhang_animal[,2])+as.numeric(overhang2_animal[,2] )
  }else{
    print("Making dataframe plant")
    mirnadf_folding<-cbind(mirnadf,foldingFigs ,"5'Cleavage"=overhang_plant[,1],"3'Cleavage"=overhang2_plant[,1] )
    output$mirnaSeqswithplots <-  DT::renderDataTable({ mirnadf_folding[,c(1,2,3,7,8,9)]},  escape = FALSE )
    overhangs_score_plant<<- as.numeric(overhang_plant[,2])+as.numeric(overhang2_plant[,2] )
  }

  foldingFigs<<-foldingFigs

  values$successStep2<-TRUE
  print(  "values$successStep2==TRUE")
  showNotification("Done. Check RNA folding Tab", type= "message")

  })## clsose button fold

##############################################
################# Expression  ################
##############################################

### select genome
##if selected
observeEvent(input$bam, {
  values$bamfilepath<<-input$bam$datapath
  print(values$bamfilepath)
})
##if uploaded
observeEvent(input$bam1, {
  values$bamfilepath<<-paste("data/bamfiles",input$bam1, sep="/")
  print(values$bamfilepath)
})
###



observeEvent(input$ButtonExp, {
  ################ Checks if previous step was succesdul
  ### Allows to run step 3 without step 2....
  shiny::validate(
    need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
    errorClass =  showNotification("Missing succesful Step 1", type= "error")
  )
  shiny::validate(
    need(values$bamfilepath!="data/bamfiles/NULL", message="missing BAM"),
    errorClass =  showNotification("Missing BAM file", type= "error")
  )


  withProgress(message = 'Calculating the expression of:', value = 0, {
    n <- 3

  #bamfilepath<-input$bam$datapath
  ###anotation files
  precsdf_dir <- input$precs$datapath
  matdf_dir <-   input$mature$datapath
  stardf_dir <- input$star$datapath
  ##




  annot<-precsdf_dir
  incProgress(1/n, detail = paste("Precursor"))

  Fcounts_prec<-featureCounts(files=values$bamfilepath, annot.ext=annot, strandSpecific=1,isPairedEnd=F,
                               allowMultiOverlap=T, countMultiMappingReads=T,
                               isGTFAnnotationFile=T,GTF.featureType="miRNA_primary_transcript",
                               GTF.attrType="ID",useMetaFeatures=F,nthreads= detectCores(all.tests = FALSE, logical = TRUE))
  annot<-matdf_dir

  incProgress(1/n, detail = paste("Mature"))
  Fcounts_mat<-featureCounts(files=values$bamfilepath, annot.ext=annot, strandSpecific=1,isPairedEnd=F,
                              allowMultiOverlap=T, countMultiMappingReads=T,
                              isGTFAnnotationFile=T,GTF.featureType="miRNA",
                              GTF.attrType="ID",useMetaFeatures=F,nthreads= detectCores(all.tests = FALSE, logical = TRUE))

   incProgress(1/n, detail = paste("Star"))
   annot<- stardf_dir
   Fcounts_star<-featureCounts(files=values$bamfilepath, annot.ext=annot, strandSpecific=1,isPairedEnd=F,
                              allowMultiOverlap=T, countMultiMappingReads=T,
                              isGTFAnnotationFile=T,GTF.featureType="miRNA",
                              GTF.attrType="ID",useMetaFeatures=F,nthreads= detectCores(all.tests = FALSE, logical = TRUE))
  # join in a single expression table
  counts_Table<-data.frame(Fcounts_prec$counts,Fcounts_mat$counts,Fcounts_star$counts )
  colnames(counts_Table)<-c("Precursor","Mature","Star")


  })  ## expression progress indicator


  Score_expression_animal<-rep(0, nrow(counts_Table)) # set all scores to zero
  Score_expression_plant<-rep(0, nrow(counts_Table)) # set all scores to zero

  ### Update scores based on expression


  for (rows in 1:nrow(counts_Table)) {
    if(counts_Table$Mature[rows] >5 ) { # more X reads matuer +1.5 points
      Score_expression_animal[rows] =  Score_expression_animal[rows]+3
      Score_expression_plant[rows] =  Score_expression_plant[rows]+4.25

    }
    if(counts_Table$Star[ rows]  >5 ) {
      Score_expression_animal[rows] = Score_expression_animal[rows]+3 # more X reads matuer +1.5 points
      Score_expression_plant[rows] = Score_expression_plant[rows]+4.25
    }
    ####### if reads precursor - reads star - reads mature is still high, more than the 10% of the mature+star, (has lots of reads on the loop!) -->  penalty!
    if(  (counts_Table$Precursor[rows] - counts_Table$Mature[rows] - counts_Table$Star[rows])  > (counts_Table$Mature[rows]+counts_Table$Star[rows])*0.8  ) {
      Score_expression_animal[rows] = Score_expression_animal[rows]-5 # more X reads matuer +1.5 points
      Score_expression_plant[rows] = Score_expression_plant[rows]-9

      print("Super penalty")
      print(counts_Table[rows,] )
    }else if((counts_Table$Precursor[rows] - counts_Table$Mature[rows] - counts_Table$Star[rows])  > (counts_Table$Mature[rows]+counts_Table$Star[rows])*0.1){
      Score_expression_animal[rows] = Score_expression_animal[rows]-5 # more X reads matuer +1.5 points
      Score_expression_plant[rows] = Score_expression_plant[rows]-8

      print("penalty")
      print(counts_Table[rows,] )
    }
  }


  #################### if input data was 5p / 3p, check whch is the mature/star
  maturesvector<-NULL
  maturesvector$Id<-rownames(counts_Table)

  if (input$matureorarm == "arm5p3p"){
      maturesvector$whoismature<-ifelse(counts_Table$Mature>counts_Table$Star, "5p is Mature", "3p is Mature" )
      maturesvector$matureis<-ifelse(counts_Table$Mature>counts_Table$Star, "5p", "3p" )
  }
  maturesvector<<-as.data.frame(maturesvector)
  ####################



  #Make it global
  Score_expression_animal<<-Score_expression_animal
  Score_expression_plant<<-Score_expression_plant

  counts_Table<<-counts_Table

  showNotification("Done Checking Expression", type= "message")
  withProgress(message = 'Loading bam...(be patient)', value = 0, {


  coverage_p<-list()
  coverage_n<-list()

    n <- 4

    incProgress(1/n, detail = "Loading bam")

    alignment<-readGAlignments(values$bamfilepath)

    incProgress(1/n, detail = "Loading bam")

    coverage_p<-coverage(alignment[strand(alignment) == "+"])

    incProgress(1/n, detail = "Loading bam")

    coverage_n<-coverage(alignment[strand(alignment) == "-"])

    incProgress(1/n, detail = "Loading bam")


  }) #close progress indicator

    withProgress(message = 'Creating Expression plots...', value = 0, {
      n<-nrow(mirnadf)

### Make it fast for trials
    #for(i in 1:3 ){
    for(i in 1:nrow(mirnadf) ){

    ##convert in GRanges
    mirname=mirnadf$ID[i]

    selectedRange <- GRanges(precsdf_adjusted$seqid[i],IRanges(start=as.numeric(precsdf_adjusted$start[i] ),end=as.numeric(precsdf_adjusted$end[i])),strand = precsdf_adjusted$strand[i])+ extrabases

    if (as.character(strand(selectedRange))=="-"){

        selectedRange_coverage<- as.numeric(unlist(coverage_n[range(selectedRange)]))

    }else{
      selectedRange_coverage<- as.numeric(unlist(coverage_p[range(selectedRange)]))
      }

   # split seq as to be used as label
    seq<-toString(mirnadf[i,]$precseqs_extended)
    seq<-unlist(strsplit(seq, split=""))


    greybars<-selectedRange_coverage
    greybars[greybars>0]<-max(selectedRange_coverage,selectedRange_coverage)
    toplot<-rbind(libs=t(selectedRange_coverage), expression=greybars)


    png(filename = paste("www/plots/",mirname, ".png",sep=""),   width = 1200, height = 480)
              par(mar=c(2,4.5,2,0))
              plot<-barplot(toplot, axes=TRUE, ylab="Number of Reads", main=mirname, col=c("blue", "grey"), border=c("blue","grey"),beside=T)
             # mtext(at = plot, text = seq,col="black", side = 1,  line = 0, cex=1)

              stardf_extrabases<-GRanges(stardf)
              matdf_extrabases<-GRanges(matdf)


              star_i=start(stardf_extrabases[i])-start(selectedRange)+1
              star_e=end(stardf_extrabases[i])-start(selectedRange)+1
              matureseq_i=start(matdf_extrabases[i])-start(selectedRange)+1
              matureseq_e=end(matdf_extrabases[i])-start(selectedRange)+1
              mtext(at = plot[1,c(1:matureseq_i,matureseq_e:length(seq))], text =seq[c(1:matureseq_i,matureseq_e:length(seq))] ,col="black", side = 1,  line = 0, cex=1)
              mtext(at = plot[1,star_i:star_e], text = seq[star_i:star_e],col="blue", side = 1,  line = 0, cex=1)
              mtext(at = plot[1,matureseq_i:matureseq_e], text = seq[matureseq_i:matureseq_e],col="red", side = 1,  line = 0, cex=1)
        dev.off()

        incProgress(1/n, detail = paste("Plot", i, "of", n))

  }#close for

    ExpressionPlot= paste("<img src=\"plots/",mirnadf$ID,".png\" width=\"1000\" height=\"600\"></img>", sep="")

    if (input$matureorarm == "maturestar"){ # if input data was mature/star
      mirnadf_plots<-data.frame("ID"=mirnadf$ID,"Mature seq"=mirnadf$mature,"Star seq"=mirnadf$star, "Reads Precursor"=counts_Table$Precursor, "Reads Mature"=counts_Table$Mature, "Reads Star"=counts_Table$Star)
      output$PLOTS <-  DT::renderDataTable({ mirnadf_plots},  escape = FALSE, selection = 'single' )
      ExpressionPlot<<-ExpressionPlot
    }else{# if was 5p 3p
      mirnadf_plots<-data.frame("ID"=mirnadf$ID,"5P arm seq"=mirnadf$mature,"3P arm seq"=mirnadf$star, "Mature"=maturesvector$matureis, "Reads Precursor"=counts_Table$Precursor, "Reads 5P arm"=counts_Table$Mature, "Reads 3P arm "=counts_Table$Star)
      output$PLOTS <-  DT::renderDataTable({ mirnadf_plots},  escape = FALSE, selection = 'single' )
      ExpressionPlot<<-ExpressionPlot
    }


     observeEvent(input$PLOTS_rows_selected, {
      row <- input$PLOTS_rows_selected
      print(row)
      output$plotouput <- renderUI({HTML(as.character(ExpressionPlot[as.numeric(row)]))})
    })

     showNotification("Done. Check Expression Plots Tab", type= "message")
     values$successStep3<-TRUE


     })# closprogress plots

})# clos button expression PLOTS


##############################################
################# Homology    ################
##############################################
observeEvent(input$Homology, {


  ### Allows to run step 4 without step 2 ,3 &4...
  shiny::validate(
    need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
    errorClass =  showNotification("Missing succesful Step 1", type= "error")
  )

  if (input$matureorarm == "arm5p3p"){
         shiny::validate(
          need( values$successStep3==TRUE, message = ('Missing succesful Step 3')),
          errorClass =  showNotification("Missing succesful Step 3", type= "error")
        )
  }else{"can proceed"}


  withProgress(message = 'Aligning matures:', value = 0, {
    n <- dim(mirnadf)[1]
  #load and prepare database (from miRBase)
  mirbase <<- scanFa(FaFile(paste("data/database","Mirbase_mature.fa", sep="/")),as="RNAStringSet")
  names(mirbase)=sapply(strsplit(as.character(names(mirbase))," "), `[`, 1)

  mirbase_organisms<-read.table("data/organisms.txt",sep="\t",header = F)## read organisms info

  if (specie=="Animal"){
    speciesforaligning<-  mirbase_organisms[grep("Metazoa",  mirbase_organisms$V4), ]
    mirNAStoalign_index<-sapply(strsplit(as.character(names(mirbase)),"-"), `[`, 1)%in%speciesforaligning$V1
    mirNAStoalign<-mirbase[mirNAStoalign_index]

  }else{
    speciesforaligning<-  mirbase_organisms[grep("Viridiplantae",  mirbase_organisms$V4), ]
    mirNAStoalign_index<-sapply(strsplit(as.character(names(mirbase)),"-"), `[`, 1)%in%speciesforaligning$V1
    mirNAStoalign<-mirbase[mirNAStoalign_index]
  }





  score_penaltylength_plant<<-NULL
  score_penaltylength_animal<<-NULL
  score2_animal<<-NULL
  score2_plant<<-NULL
  conservationtype<<-NULL

  Alignments<-function(datadf){
    n=nrow(mirnadf)
    incProgress(1/n, detail = paste("Alignments"))
    toreturn<-NULL
    #print(datadf[7])

    if (input$matureorarm == "maturestar"){  # if input data has  mature/arm info
        seed1<-as.character(RNAStringSet(DNAStringSet(substr(as.character(datadf[2]), 2,9))))
        matureasRNA<<-as.character(RNAStringSet(DNAStringSet(as.character(datadf[2]))))
    }else{# if I computed mature/arm using expression data
      #head(mirnadf)
      matureasRNA_0<<-ifelse(maturesvector[maturesvector$Id==datadf[1], ]$matureis =="5p", datadf[2] , datadf[3])
      seed1<-as.character(RNAStringSet(DNAStringSet(substr(as.character(matureasRNA_0), 2,9))))
      matureasRNA<<-as.character(RNAStringSet(DNAStringSet(as.character(matureasRNA_0))))
    }



    toalignwhole<-mirNAStoalign[grep(matureasRNA, as.character(mirNAStoalign)), ] ### select identical mirnas
    toalignseed<-mirNAStoalign[grep(paste("^.",seed1,sep=""), as.character(mirNAStoalign)), ] ### select identical mirnas

    #names(matureasRNA)<-paste0("<div> <span style=\"color:red\"><strong>CANDIDATE </strong> </span></div>")
    names(matureasRNA)<-"CANDIDATE"

    toalign<-RNAStringSet(matureasRNA)

    numberofalignments<-20

    if (length(toalignwhole)>1 ){ # if some are dientical hits
      if (length(toalignwhole)<numberofalignments &  (length(toalignseed)+length(toalignwhole))>numberofalignments  ){ # if identical are less than numberofalignments but plus same seed more than numberofalignments
          toalign<-c(toalign,toalignwhole )####
          toalign<-c(toalign, sample(toalignseed,numberofalignments-length(toalignwhole)) )#### we ad same seed until numberofalignments

          if(length(toalignwhole)>6){## if there are between 5 and 20 identicals  (+ >20 similar)
            score2_animal<<-c(score2_animal,3.5)
            score2_plant<<-c(score2_plant,6)
            conservationtype<<-c(conservationtype,"strong 6")
          }else{
            score2_animal<<-c(score2_animal,3)###### if there are between 2-4 identicals (+ >20 similar)
            score2_plant<<-c(score2_plant,4)
            conservationtype<<-c(conservationtype,"strong 4")

          }
      }

      if (length(toalignwhole)<numberofalignments & length(toalignwhole)>0 & length(toalignseed)+length(toalignwhole)<numberofalignments  ){
        toalign<-c(toalign, toalignwhole)####if identical + same seed less numberofalignments, we use all
        toalign<-c(toalign, toalignseed)


         if(length(toalignwhole)>5){## if there are between 5 and 20 identicals  (+ <20 similar)
          score2_animal<<-c(score2_animal,3.4)
          score2_plant<<-c(score2_plant,5.8)
          conservationtype<<-c(conservationtype,"strong 5.8")

        }else{
          score2_animal<<-c(score2_animal,2.78)###### if there are between 2-4 identicals (+ <20 similar)
          score2_plant<<-c(score2_plant,4.8)
          conservationtype<<-c(conservationtype,"medium 4.8")

        }

      }
      if(length(toalignwhole)>numberofalignments){#### if more than numberofalignments identicals, we select ALL identical
        toalign<-c(toalign,toalignwhole)
        score2_animal<<-c(score2_animal,4)######
        score2_plant<<-c(score2_plant,6.5)######
        conservationtype<<-c(conservationtype,"very strong 6.5")


      }

    }else if(length(toalignseed)>1){#### if no identicals but some with same seed
        if (length(toalignseed)<numberofalignments  ){#if less than numberofalignments we use all
          toalign<-c(toalign,toalignseed)
          score2_animal<<-c(score2_animal,2)######
          score2_plant<<-c(score2_plant, 0.5)######
          conservationtype<<-c(conservationtype,"low 2")


        }
        if (length(toalignseed)>numberofalignments  ){#if more same seed than numberofalignments, we select ALL
          toalign<- c(toalign,toalignseed)###
          score2_animal<<-c(score2_animal,2.5)######
          score2_plant<<-c(score2_plant,0.8)######
          conservationtype<<-c(conservationtype,"low 0.8")


        }

    }

  if (length(toalign)>1){
     toalign<-toalign[!duplicated(names(toalign)),]
     Alignmentout <- msa(toalign, "ClustalOmega",order= "input")  ## Align the miRNAs
     toreturn<- c(toreturn, Alignmentout)
    }else{ toreturn<-c(toreturn,"No")
        score2_animal<<-c(score2_animal,0)
        score2_plant<<-c(score2_plant,-2.5)
        conservationtype<<-c(conservationtype,"none -2.5")


              }
        return(toreturn)
    }


  Alignmentlist<<-unlist(apply(mirnadf, 1, Alignments))

  ### Sum the score of the homology to the global score
  Score_homology_animal<<-score2_animal
  Score_homology_plant<<-score2_plant


  if (input$matureorarm == "maturestar"){  # if input data has  mature/arm info
    toprint<-mirnadf[,-c(4,6)  ]
    colnames(toprint)[2]<-"Mature"
    colnames(toprint)[3]<-"Star"
    toprintaligns<<-toprint

  }else{

      toprint<-cbind(mirnadf[,-c(4,6) ],maturesvector[,3] )
      colnames(toprint)[2]<-"5p arm"
      colnames(toprint)[3]<-"3p arm"
      colnames(toprint)[5]<-"Mature"
      toprintaligns<<-toprint

  }

  })#close progress bar


  showNotification("Done. Check Homology Tab", type= "message")
  values$successStep4<-TRUE


  #create reactive output
  print(head(toprintaligns))
  output$alignmentsoutput<-DT::renderDataTable({ toprintaligns },selection = "single",server = FALSE,
                                              options = list(paging=TRUE,
                                                             searching=FALSE,
                                                             filtering=FALSE,
                                                             ordering=TRUE))


  observeEvent(input$alignmentsoutput_rows_selected, {
    row <- input$alignmentsoutput_rows_selected
    if(class(Alignmentlist[[as.numeric(row)]])=="character"){## if has no alignments
      output$alignments2 <- renderPrint(Alignmentlist[[as.numeric(row)]] )
    }else{#if there are lignments
      alignemntsasdf<-as.data.frame(Alignmentlist[[as.numeric(row)]]@unmasked)
      colnames(alignemntsasdf)<-c("")
      output$alignments2 <- renderPrint(alignemntsasdf)
    }
  })# close reactive output


  })#close button homology


##############################################
################# Integration ################
##############################################
observeEvent(input$ButtonIntegrate, {

  ### Must hav ALL steps to be able to  calculate the real score
  shiny::validate(
    need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
    errorClass =  showNotification("Missing succesful Step 1", type= "error")
  )
  shiny::validate(
    need( values$successStep2==TRUE, message = ('Missing succesful Step 2')),
    errorClass =  showNotification("Missing succesful Step 2", type= "error")
  )
  shiny::validate(
    need( values$successStep3==TRUE, message = ('Missing succesful Step 3')),
    errorClass =  showNotification("Missing succesful Step 3", type= "error")
  )
  shiny::validate(
    need( values$successStep4==TRUE, message = ('Missing succesful Step 4')),
    errorClass =  showNotification("Missing succesful Step 4", type= "error")
  )



  if(specie=="Animal"){
    Score<-Score_expression_animal+overhangs_score_animal+Score_homology_animal+penaltyscorelength_animal

  }else{
    Score<-Score_expression_plant+overhangs_score_plant+Score_homology_plant+penaltyscorelength_plant

  }

   mirnadf_integrated<-cbind(mirnadf,counts_Table,conservationtype,Score)
   mirnadf_integrated_globalcopie<<-cbind(mirnadf,counts_Table,conservationtype,Score)

  print("Integrating")


 # mirnadf_integrated$box<-shinyInput(checkboxInput, nrow(mirnadf_integrated), 'v2_', value = TRUE)
  shinyInput = function(FUN, len, id, ...) {
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...))
    }
    inputs
  }

  res = data.frame(
    Name=mirnadf_integrated[,1],
    NewName = shinyInput(textInput, nrow(mirnadf_integrated), 'v1_', value = ""),
    Filter = shinyInput(checkboxInput, nrow(mirnadf_integrated), 'v2_', value = FALSE,width=0.5),# Checkbox: True=check False= not check
    mirnadf_integrated[,c(-1,-5,-6)],#[,c(2,3,4,7,8,9,10)],
    stringsAsFactors = FALSE
  )


  ##### SET score threshold!
  scorethreshold<<-14
  #####
## pre-check some rows absed on their score
  res$NewName<- paste0('<div class=\"form-group shiny-input-container\">\n  <input id=\"v1_1\" type=\"text\" class=\"form-control\" value=\"', mirnadf_integrated$ID,'">\n</div>')
  res[res$Score>=scorethreshold,]$Filter<- gsub('\"checkbox\"/>\n', '\"checkbox\" checked=\"checked\"/>\n',  res[res$Score>=scorethreshold,]$Filter )
####
  res_globalcopie<<-res


  if (input$matureorarm == "maturestar"){ ## if input data has mature/star info
      output$Intergartiontable = DT::renderDataTable(
        res, server = FALSE, escape = FALSE, selection = 'single', rownames= FALSE, options = list(
          preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
          drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } '),
          pageLength = 5
        )
      )
  }else{## if mature/star was computed in mirplot

    res_v2<-cbind(res[,1:5], maturesvector[,3], res[ , -1:-5]  )
    colnames(res_v2)[4]<-"5p Arm"
    colnames(res_v2)[5]<-"3p Arm"
    colnames(res_v2)[6]<-"Mature"

    colnames(res_v2)[9]<-"Reads 5p arm"
    colnames(res_v2)[10]<-"Reads 3p arm"

    output$Intergartiontable = DT::renderDataTable(
      res_v2, server = FALSE, escape = FALSE, selection = 'single', rownames= FALSE, options = list(
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } '),
        pageLength = 5
      )
    )

  }

  # print the values of inputs
  shinyValue = function(id, len) {
    unlist(lapply(seq_len(len), function(i) {
      value = input[[paste0(id, i)]]
      if (is.null(value)) NA else value
    })) }


###############
### react on selecting rows (Not in checkbox)
  observeEvent(input$Intergartiontable_rows_selected, {
    rowInt <- input$Intergartiontable_rows_selected
    print(rowInt)
    output$Intergartiontable2 <- renderUI({HTML(as.character(foldingFigs[as.numeric(rowInt)]))})
    output$Intergartiontable3 <- renderUI({HTML(as.character(ExpressionPlot[as.numeric(rowInt)]))})
    #output$Intergartiontable4 <- renderPrint(Alignmentlist[as.numeric(rowInt)])


    if(class(Alignmentlist[[as.numeric(rowInt)]])=="character"){## if has no alignments
      output$Intergartiontable4 <- renderPrint(Alignmentlist[[as.numeric(rowInt)]] )
    }else{#if there are lignments
      alignemntsasdf<-as.data.frame(Alignmentlist[[as.numeric(rowInt)]]@unmasked)
      colnames(alignemntsasdf)<-c("")
      output$Intergartiontable4 <- renderPrint(alignemntsasdf)
    }



  })

  rowSelect <- reactive({
    paste(sort(unique(input[["rows"]])),sep=',')
    print(paste("Checked", rowSelect)    )
  })




###############
############### Download data
  # Downloadable csv of SELECTED ROWS ----
  output$downloadData <- downloadHandler(
    filename = "FilteredMiRNAs.csv",
    content = function(file) {

      ## first, get the pre-select  rows
      originalstatus<- data.frame(v1 = res_globalcopie$Name ,
                                   v2 = res_globalcopie$Score )
      originalstatus$v2<-ifelse(originalstatus$v2>=scorethreshold, TRUE, FALSE )
      originalstatus<<-originalstatus
      #T2nd we get the info o fthe items that user modified
      Selectiondatframe<<- data.frame(v1 = shinyValue('v1_', nrow(mirnadf_integrated)),
                                     v2 = shinyValue('v2_', nrow(mirnadf_integrated)))
      Selectiondatframe$v1<-as.character(Selectiondatframe$v1)

      # if not modified, has NA, and therfore we assign it the value of the orifinal pre-selection
      for(i in 1:nrow(Selectiondatframe)){
        if(is.na(Selectiondatframe[i,1])){
          Selectiondatframe[i,1]<-as.character(originalstatus[i,1])
        }
        if(is.na(Selectiondatframe[i,2])){
          Selectiondatframe[i,2]<-originalstatus[i,2]
        }
      }

      Todownload<<-mirnadf_integrated[Selectiondatframe$v2==TRUE,1:5]
      Todownload<<-cbind(Selectiondatframe[Selectiondatframe$v2==TRUE,1],Todownload)
      write.csv(Todownload, file, row.names = FALSE)
    }
  )



  # Downloadable MATURES of SELECTED ROWS ----
  output$downloadMatures <- downloadHandler(
    filename = "Maturesafterfilter.fa",
         content = function(file) {
           if (input$matureorarm == "maturestar"){ ## if input data has mature/star info

              ## first, get the pre-select  rows
              originalstatus<<- data.frame(v1 = res_globalcopie$Name ,
                                           v2 = res_globalcopie$Score )
              originalstatus$v2<-ifelse(originalstatus$v2>=scorethreshold, TRUE, FALSE )

              #T2nd we get the info o fthe items that user modified
              Selectiondatframe<<- data.frame(v1 = shinyValue('v1_', nrow(mirnadf_integrated)),
                                              v2 = shinyValue('v2_', nrow(mirnadf_integrated)))
              Selectiondatframe$v1<-as.character(Selectiondatframe$v1)

              # if not modified, has NA, and therfore we assign it the value of the orifinal pre-selection
              for(i in 1:nrow(Selectiondatframe)){
                if(is.na(Selectiondatframe[i,1])){
                  Selectiondatframe[i,1]<-as.character(originalstatus[i,1])

                }
                if(is.na(Selectiondatframe[i,2])){
                  Selectiondatframe[i,2]<-originalstatus[i,2]
                }
              }

              Todownload<<-mirnadf_integrated[Selectiondatframe$v2==TRUE,1:5]
              Todownload<<-cbind(Selectiondatframe[Selectiondatframe$v2==TRUE,1],Todownload)
              Todownloadasfasta<-paste0(">",Todownload[,1],"\n",Todownload[,3])
              write(Todownloadasfasta, file )

               }else{  ### if input was 5p / 3p

            ## first, get the pre-select  rows
              originalstatus<<- data.frame(v1 = res_globalcopie$Name ,
                                           v2 = res_globalcopie$Score )
              originalstatus$v2<-ifelse(originalstatus$v2>=scorethreshold, TRUE, FALSE )

              #T2nd we get the info o fthe items that user modified
              Selectiondatframe<<- data.frame(v1 = shinyValue('v1_', nrow(mirnadf_integrated)),
                                              v2 = shinyValue('v2_', nrow(mirnadf_integrated)))
              Selectiondatframe$v1<-as.character(Selectiondatframe$v1)

              # if not modified, has NA, and therfore we assign it the value of the orifinal pre-selection
              for(i in 1:nrow(Selectiondatframe)){
                if(is.na(Selectiondatframe[i,1])){
                  Selectiondatframe[i,1]<-as.character(originalstatus[i,1])

                }
                if(is.na(Selectiondatframe[i,2])){
                  Selectiondatframe[i,2]<-originalstatus[i,2]
                }
              }


              ## maturesvector$matureis
              ### check if mature is 5p or 3p and return mature!
              Todownload_0<<-cbind(mirnadf_integrated,Selectiondatframe, matureis=maturesvector$matureis )
              Todownload_1<<-Todownload_0[Todownload_0$v2==TRUE,c(1:5,13 )]

              Todownloadasfasta2<<- ifelse(Todownload_1$matureis=="5p", paste0(">",Todownload_1[,1],"\n",Todownload_1[,2]),  paste0(">",Todownload_1[,1],"\n",Todownload_1[,3]))

              #Todownloadasfasta<-paste0(">",Todownload[,1],"\n",Todownload[,3])
              write(Todownloadasfasta2, file )

              }
    }#content end
 )##end download handler
  toreturnall<<-NULL
  if (input$matureorarm == "maturestar"){
    toreturnall<<-mirnadf_integrated_globalcopie
    #write.csv(toreturnall, file3 )
  }else{
    toreturnall<<-cbind(mirnadf_integrated_globalcopie, maturesvector)
    colnames(toreturnall)[2]<-"5p arm"
    colnames(toreturnall)[3]<-"3p arm"
    colnames(toreturnall)[8]<-"5p reads"
    colnames(toreturnall)[9]<-"3p reads"
    #write.csv(toreturnall, file3 )
  }



  # Download ALL  ----
  output$downloadALL <- downloadHandler(
    filename = "All_data.txt",
    content = function(file3) {
      write.csv(toreturnall, file3 )
    }#content end
  )##end download ALL handler


  ##################### Score stats

  # SCORES Histogram overlaid with kernel density curve
  output$Scoreshistogram<- renderPlot({ ggplot(toreturnall, aes(x=Score)) +
                      geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                      binwidth=.5,
                      colour="black", fill="white") +
                    geom_density(alpha=.2, fill="#FF6666")+theme(legend.position="none",axis.text=element_text(size=13), axis.title=element_text(size=15)) }) # Overlay with transparent density plot

   # SCORES boxplot/cateogry
  Correct<<-NULL
  Doubtful<<-NULL
  False<<-NULL
    ## I need to put If to avoid crash if 0 rows selected...
  if( nrow(toreturnall[toreturnall$Score>=scorethreshold,])>0 ){
    Correct<<-data.frame(toreturnall[toreturnall$Score>=scorethreshold,c(7,8,9,10,11)],"class"=paste('Correct miRNAs (Scores [20,',scorethreshold,'])'))
    colnames(Correct)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }else{
    Correct<<-data.frame(0,0,0,"A",0,paste('Correct miRNAs (Score= [20,',scorethreshold,'])'))
    colnames(Correct)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }
  if( nrow(toreturnall[toreturnall$Score<scorethreshold & toreturnall$Score>=scorethreshold-6 ,])>0 ){
    Doubtful<<-data.frame(toreturnall[toreturnall$Score<scorethreshold & toreturnall$Score>=scorethreshold-6 ,c(7,8,9,10,11)],"class"=paste('Doubtful miRNAs (Scores (',scorethreshold,',',scorethreshold-6,'])'))
    colnames(Doubtful)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }else{
    Doubtful<<-data.frame(0,0,0,"A",0,paste('Doubtful miRNAs (Scores (',scorethreshold,',',scorethreshold-6,'])'))
    colnames(Doubtful)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }
  if( nrow(toreturnall[toreturnall$Score<scorethreshold-6 ,])>0 ){
    False<<-data.frame(toreturnall[toreturnall$Score<scorethreshold-6 ,c(7,8,9,10,11)],"class"=paste('False miRNAs (Scores (',scorethreshold-6,',-10])'))
    colnames(False)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }else{
    False<<-data.frame(0,0,0,"A",0,paste('False miRNAs (Scores (',scorethreshold-6,',-10])'))
    colnames(False)<-c(colnames(toreturnall[,c(7,8,9,10,11)]),"class")
  }
  Scoresboxplotdata <<-NULL
  Scoresboxplotdata <<-base::rbind(Correct,Doubtful,False)
  colnames(Scoresboxplotdata)<-colnames(Correct)

  output$ScoreBoxplots1<- renderPlot({ggplot(Scoresboxplotdata, aes(x = class , y = Score)) +
    stat_summary(fun.y=mean, colour="red", geom="point") +
    geom_boxplot (aes(fill=Score), alpha=.5, width=1, position = position_dodge(width = 1),  outlier.colour = "dark gray", outlier.size = 1)+
    ggtitle("Scores per class") + theme(legend.position="none",axis.text=element_text(size=13), axis.title=element_text(size=15))})

  output$ScoreBoxplots2<- renderPlot({

    ggplot(Scoresboxplotdata, aes(x = class , y = log(Score))) +
      stat_summary(fun.y=mean, colour="red", geom="point") +
      geom_boxplot (aes(fill=Score), alpha=.5, width=1, position = position_dodge(width = 1),  outlier.colour = "dark gray", outlier.size = 1)+
      ggtitle("Log(Scores per class)") + theme(legend.position="none",axis.text=element_text(size=13), axis.title=element_text(size=15))
    })




  })# close button integration

}#close server
