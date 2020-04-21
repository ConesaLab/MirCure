# Define server logic to summarize and view selected dataset ----

###### In cluster, activate viennarna module.
#system("module load viennarna")

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
library("rmarkdown")

memory.size(max=TRUE)

options(shiny.maxRequestSize = 10000*1024^2)


##############################################
############ Score definition ################
##############################################
## Overhang scores
perfectmatch_animal<<-1
acceptable_animal<<-0.75
suspicious_animal<<-0.5

perfectmatch_plant<<-1
acceptable_plant<<-0
suspicious_plant<<-0


Nolooppenalty<<-(-5) # if there is no loop at all, is bad for sure
# bad & uncommon structures  = 0


## miRNA arm lengths
penalty_length<<-(-2) # penalty applied for  3p and 5p if <18 or >23

## Expression thresholds (applied to each arm)
#more reads than ExpLevel1
ExpLevel1<<- 5 #minim reads
ExpLevel1_pointsAnimal<<-1.5
ExpLevel1_pointsPlant<<-5

#more reads than ExpLevel2
ExpLevel2<<-100 #minim reads
ExpLevel2_pointsAnimal<<-2
ExpLevel2_pointsPlant<<-6


#Too many reads in loop has a penalty
loopreadthreshold<<-0.1  #if loop contains more than X% of the pre-miRNA reads
penatlytoomanyreadsloopAnimal<<-(5)
penatlytoomanyreadsloopPlant<<-(5)

#Too many reads in loop has a penalty
### Conservation scores (important for Plants!)

#### >20 miRNAs in mirbase identical
Score_conservation_20id_Animal<<-1
Score_conservation_20id_Plant<<-2

#### 6-20 miRNAs identicals, + similar (same seed) >20
Score_conservation_6_20id20_Animal<<-1
Score_conservation_6_20id20_Plant<<-2

#### 2-5 miRNAs identicals, + similar (same seed) >20
Score_conservation_2_5id20_Animal<<-0.75
Score_conservation_2_5id20_Plant<<-1.5

#### 6-20 miRNAs identicals, + similar (same seed) <20
Score_conservation_6_20id_Animal<<-0.5
Score_conservation_6_20id_Plant<<-1

#### 2-5 miRNAs identicals, + similar (same seed) <20
Score_conservation_2_5id_Animal<<- (0.25)
Score_conservation_2_5id_Plant<<-0.5

#### 0 identicals, + similar (same seed) >20
Score_conservation_0id20_Animal<<- (0)
Score_conservation_0id20_Plant<<- 0.5


#### 0 identicals, + similar (same seed) <20
Score_conservation_0id_Animal<<-(-0.25)
Score_conservation_0id_Plant<<- -0.5


#### Final Score Formula
CalculateScore<<-function(expression, overhang,homology,penaltylen){
  finalscore<<-expression+overhang+homology+penaltylen
  return(finalscore)
}


folded_globe <<- list()

###########################################


server <- function(input, output, session) {
  #######
  ###### variables that are session-specific
  #vars to track progress
  values <- reactiveValues()

  values$successStep1<-FALSE
  values$successStep2<-FALSE
  values$successStepAdjust<-FALSE
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
    #precsdf <- read.table("")


    return(head(precsdf))
  })

  #mat file
  output$mature <- renderTable({
    req(input$mature)
    matdf <- readGFF(input$mature$datapath)
    #matdf <- readGFF("")


    return(head(matdf))
  })

  #star file
  output$star <- renderTable({
    req(input$star)
    stardf <- readGFF(input$star$datapath)
    #stardf <- readGFF(")
    return(head(stardf))
  })

  ## user uplaoded mature/star or 5p/3p

  observeEvent(input$matureorarm, {
    if(input$matureorarm == "maturestar"){
      maturestar=TRUE
      color_arm1<<-"red"
      color_arm2<<-"blue"
      print("maturestar")}

    if(input$matureorarm == "arm5p3p"){
      maturestar=FALSE
      color_arm1<<-"darkorange3"
      color_arm2<<-"darkgreen"
      print("arm5p3p")}
  })

  #################################
  ######### Get sequences #########
  #################################

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


    #genomefasta <- FaFile(paste("data/genomes","genome.fasta", sep="/"))
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

    matseqs <<- getSeq(genomefasta,GRanges(matdf$seqid,IRanges(start=as.numeric(matdf$start),end=as.numeric(matdf$end)),strand =precsdf$strand ))
    starseqs <<- getSeq(genomefasta,GRanges(stardf$seqid,IRanges(start=as.numeric(stardf$start),end=as.numeric(stardf$end)),strand =precsdf$strand ))
    ##########

    #extrabases<<-11
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
    precseqs <<- getSeq(genomefasta,GRanges(precsdf_adjusted$seqid,IRanges(start=as.numeric(precsdf_adjusted$start ),end=as.numeric(precsdf_adjusted$end)),strand =precsdf_adjusted$strand ) )



    #output$mirnaSeqs <-renderTable({
    #mirnadf<<-data.frame(ID=as.character(precsdf$ID), mature=as.character(matseqs),star=as.character(starseqs),"Loci"=paste(precsdf$seqid,":",precsdf$start,"-",precsdf$end,":",precsdf$strand, sep=''),precursor=as.character(precseqs),precseqs_extended=as.character(precseqs_extended) )
    colnames(precsdf) = c("seqid" ,"source", "type",   "start",  "end",    "score"  ,"strand", "phase"  ,"ID" )
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

    ########### penalty nature and star lengths

    ########


  })## clsose button seqs
  #################################
  ##### Check arms annotation #####
  #################################

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



  observeEvent(input$ButtonCheck, {
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



    ####################

    showNotification("Checking Expression", type= "message")
    withProgress(message = 'Loading bam...(be patient)', value = 0, {

      n <- 4
      coverage_p<<-list()
      coverage_n<<-list()


      incProgress(1/n, detail = "Loading bam")

      #alignment <<- readGAlignments(".bam")

       alignment <<- readGAlignments(values$bamfilepath)

      incProgress(1/n, detail = "Loading bam")

      coverage_p<<-coverage(alignment[strand(alignment) == "+"])

      incProgress(1/n, detail = "Loading bam")

      coverage_n<<-coverage(alignment[strand(alignment) == "-"])

      incProgress(1/n, detail = "Loading bam")


    }) #close progress indicator

    withProgress(message = 'Adjust structure', value = 0, {


      ### Make it fast for trials

      adjustStarPosition <<-  list()
      adjustMaturePosition <<- list()
      adjustMatureSequence <<- c()
      adjustStarSequence <<- c()
      n <- 0
      q <- nrow(mirnadf)

      for(i in 1: nrow(mirnadf) ){
        incProgress(1/q, detail = paste("Adjusting the annotation:", i, "of", q))

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


        stardf_extrabases<-GRanges(stardf)
        matdf_extrabases<-GRanges(matdf)


        star_i <- start(stardf_extrabases[i])-start(selectedRange)+1
        star_e <- end(stardf_extrabases[i])-start(selectedRange)+1
        matureseq_i <- start(matdf_extrabases[i])-start(selectedRange)+1
        matureseq_e <- end(matdf_extrabases[i])-start(selectedRange)+1



        ##############correct  annotation based on the expression


        starReads <- mean(selectedRange_coverage[star_i:star_e]) # average reads in the previous annotation
        matureReads <- mean(selectedRange_coverage[matureseq_i:matureseq_e]) # average reads in the previous annotation

        ## There is no reads, we just use the previous annotation
        if (starReads == 0 | matureReads == 0) {
          adjustStarPosition[[i]] <- c(star_i, star_e) # get the position (star/end) for check structure
          adjustMaturePosition[[i]] <- c (matureseq_i, matureseq_e)
          adjustMatureSequence[i] <- as.character(mirnadf$mature[i])
          adjustStarSequence[i] <- as.character(mirnadf$star[i])
        } else {
          readData <- as.data.frame(cbind(seq, as.numeric(selectedRange_coverage)))

          TrueStar <- c()
          TrueMature <- c()

          colnames(readData) <- c ("nucleotide", "position")
          for (j in (star_i -4) : (star_e + 4)) {
            # We only go into this if-else if we could find [starReads * 0.8< reads <starReads * 1.6]
            if ( as.numeric(as.character(readData$position[j])) > starReads * 0.8 &
                 as.numeric(as.character(readData$position[j])) < starReads * 1.6 ){
              if (j == star_i-4) {
                continueReads1 <- star_i - 4
              } else if ( j == n + 1) {
                continueReads1 <- append(continueReads1,j)
              } else {
                continueReads1 <- j
              }
              if (length(continueReads1) >= length(TrueStar)) {
                TrueStar <- continueReads1
              }
              n <- j # save last number
            }
          }


          for (k in (matureseq_i - 4) : (matureseq_e + 4)) { ## We
            if ( as.numeric(as.character(readData$position[k])) > matureReads * 0.8 &
                 as.numeric(as.character(readData$position[k])) < matureReads * 1.6 ){
              if (k == matureseq_i-4) {
                continueReads2 <- matureseq_i - 4
              } else if ( k == n + 1) {
                continueReads2 <- append(continueReads2,k)
              } else {
                continueReads2 <- k
              }
              if (length(continueReads2) >= length(TrueMature)) {
                TrueMature <- continueReads2
              }
              n <- k # save last number
            }
          }



          ## We should ensure the adjust sequence in the range of 20 to 26
          if (length(TrueStar) >= 20 && length(TrueStar)<= 26 && length(TrueMature) >= 20 && length(TrueStar) <= 26) {
            adjustStarPosition[[i]] <- c(min(TrueStar), max(TrueStar)) # get the position (star/end) for check structure
            adjustMaturePosition[[i]] <- c (min(TrueMature), max(TrueMature))
            adjustMatureSequence[i] <- paste(as.character(readData$nucleotide[min(TrueMature):max(TrueMature)]), collapse = '')
            adjustStarSequence[i] <- paste(as.character(readData$nucleotide[min(TrueStar):max(TrueStar)]), collapse = '')
          } else if ((length(TrueStar) >= 20 && length(TrueStar) <= 26) && (length(TrueMature) < 20 | length(TrueStar) >26 )){
            adjustStarPosition[[i]] <- c(min(TrueStar), max(TrueStar)) # get the position (star/end) for check structure
            adjustMaturePosition[[i]] <- c (matureseq_i, matureseq_e)
            adjustMatureSequence[i] <- as.character(mirnadf$mature[i])
            adjustStarSequence[i] <- paste(as.character(readData$nucleotide[min(TrueStar):max(TrueStar)]), collapse = '')
          } else if ((length(TrueMature) >= 20 && length(TrueMature) <= 26) && (length(TrueStar) <20 | length(TrueStar) > 26)) {
            adjustStarPosition[[i]] <- c(star_i, star_e) # get the position (star/end) for check structure
            adjustMaturePosition[[i]] <- c (min(TrueMature), max(TrueMature))
            adjustMatureSequence[i] <- paste(as.character(readData$nucleotide[min(TrueMature):max(TrueMature)]), collapse = '')
            adjustStarSequence[i] <- as.character(mirnadf$star[i])
          } else {
            adjustStarPosition[[i]] <- c(star_i, star_e) # get the position (star/end) for check structure
            adjustMaturePosition[[i]] <- c(matureseq_i, matureseq_e)
            adjustMatureSequence[i] <- as.character(mirnadf$mature[i])
            adjustStarSequence[i] <- as.character(mirnadf$star[i])
          }

          continueReads1 <- c()
          continueReads2 <- c()
        }
        ####################################################

      }#close for
      adjustStarPosition <<- adjustStarPosition  # get the position (star/end) for check structure
      adjustMaturePosition <<- adjustMaturePosition
      adjustMatureSequence <<- adjustMatureSequence
      adjustStarSequence <<- adjustStarSequence


    })
    values$successStepAdjust<-TRUE
  })# clos button expression PLOTS

                                                                    #
  #################################
  ### Check secondary structure ###
  #################################

  observeEvent(input$ButtonFold, {

    ################ Checks if previous step was successful
    shiny::validate(
      need( values$successStep1==TRUE, message = ('Missing successful Step 1')),
      errorClass =  showNotification("Missing successful Step 1", type= "error")
    )
    shiny::validate(
      need( values$successStepAdjust==TRUE, message = ('Missing successful Step 2')),
      errorClass =  showNotification("Missing successful Step 2", type= "error")
    )

    #### requires Vienna, downalaod from https://www.tbi.univie.ac.at/RNA/index.html#download
    ########## ./configure  make , sudo make install
    withProgress(message = 'Folding precursors', value = 0, {
      overhang_animal<-NULL
      overhang_animal_adjust <- NULL
      overhang2_animal<-NULL
      overhang2_animal_adjust <- NULL
      overhang_plant<-NULL
      overhang_plant_adjust <- NULL
      overhang2_plant<-NULL
      overhang2_plant_adjust <-NULL
      mirnadf_mismatch<-NULL
      precsdf_1 <- as.data.frame(precsdf)
      finalStarPosition <- list()
      finalMaturePosition <- list()
      q <- nrow(mirnadf)
      penaltyscore1 <- c()
      penaltyscore2 <-c()
      penaltyscorelength_animal <-c()




      for (i in 1:nrow(mirnadf)){


        incProgress(1/q, detail = paste0("May take a long time ... ", i, " of ", q))
        #if (!str_count(mirnadf$precursor[i], "N") > 20){## if miRNA precursorcontains > 20 Ns


        #folded=run_RNAfold(as.character(mirnadf$precseqs_extended[i]), RNAfold.path = "RNAfold", detectCores(all.tests = FALSE, logical = TRUE))
        folded=run_RNAfold(as.character(mirnadf$precseqs_extended[i]), RNAfold.path = "RNAfold", parallel.cores= ifelse(detectCores(all.tests = FALSE, logical = TRUE)-1 <1,1, detectCores(all.tests = FALSE, logical = TRUE)-1))#Use all available cores minus 1, if only 1, use 1.

        coord=ct2coord(makeCt(folded[2,], folded[1,]))

        ################# Count Overhang within FOLDING ############################################
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
        maturecord1<-range(maturecord0,maturecord0+nchar(as.character(mirnadf$mature[i]))-1)

        ## we would get different position if we start from +/-
        #if(precsdf_1$strand[i] == "+") {
        maturecord1_adjust <- adjustMaturePosition[[i]]
        starcord1_adjust <- adjustStarPosition[[i]]
        #} else {
        #  MatureEnd <- str_length(mirnadf$precseqs_extended[i]) - adjustMaturePosition[[i]][1] + 1
        #  MatureStart <- str_length(mirnadf$precseqs_extended[i]) - adjustMaturePosition[[i]][2] + 1
        #  maturecord1_adjust <- c(MatureStart, MatureEnd)
        #  adjustMaturePosition[[i]] <- c(MatureStart, MatureEnd)

        #  starEnd <- str_length(mirnadf$precseqs_extended[i]) - adjustStarPosition[[i]][1] + 1
        #  starStart <- str_length(mirnadf$precseqs_extended[i]) - adjustStarPosition[[i]][2] + 1
        #  starcord1_adjust <- c(starStart, starEnd)
        #  adjustStarPosition[[i]] <- c(starStart, starEnd)
        #}
        if(maturecord1_adjust[1] > starcord1_adjust[1]) {
          temp <- maturecord1_adjust
          maturecord1_adjust <- starcord1_adjust
          starcord1_adjust <- temp
        }

        if (maturecord0[1] > starcord0[1]) {
          temp <- maturecord1
          maturecord1 <- starcord1
          starcord1 <- temp
        }



        ###################################################
        if (!str_count(mirnadf$precursor[i], "N") > 20){## if miRNA precursorcontains > 20 Nt
          ##if 5P mature
          if(maturecord1[1]<starcord1[1]){# if 5p
            print(paste("mature is 5'",i))

            overlap=maturecord1[2]-starcord1[1]+1
            if(overlap<0) {#if mature and star don't overlap (as it should)
              colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1)) ), rep(color_arm1,length(seq(maturecord1[1], maturecord1[2]))),rep("Black",length(seq(maturecord1[2]+1, starcord1[1]-1))),rep(color_arm2,length(seq(starcord1[1], starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
              foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

              foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
              foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
            }else{# If mature and star overlap (they shouldn't!) don't crash do :
              overlaFLAG=TRUE
              if (overlap>0){# if there is overlap
                colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1))), rep(color_arm1,length(seq(maturecord1[1], maturecord1[2]-overlap))),rep("Orange",length(seq(starcord1[1], maturecord1[2]))),rep(color_arm2,length(seq(starcord1[1]+overlap, starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
              if (overlap==0){# if they dont overlap, but there is no loop

                colorvector<-c(rep("Black", length(seq(1,maturecord1[1]-1)) ), rep(color_arm1,length(seq(maturecord1[1], maturecord1[2]-overlap))),rep(color_arm2,length(seq(starcord1[1]+overlap, starcord1[2]))), rep("Black",length(seq(starcord1[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
            }



          }else{   ##################### Mature 3'
            print(paste("mature is 3'", i))

            overlap=starcord1[2]+1-maturecord1[1]

            if(overlap<0){ # if no overlap

              colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep(color_arm2,length(seq(starcord1[1], starcord1[2]))),rep("Black",length(seq(starcord1[2]+1, maturecord1[1]-1))),rep(color_arm1,length(seq(maturecord1[1], maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
              foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

              foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
              foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)

            }else{# If mature and star overlap (they shouldn't!) don't crash do :

              if (overlap>0){
                colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep(color_arm2,length(seq(starcord1[1], starcord1[2]-overlap))),rep("Orange",length(seq(maturecord1[1], starcord1[2]))),rep(color_arm1,length(seq(maturecord1[1]+overlap, maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
              if (overlap==0){# if overlap is zero but also no gap
                colorvector<-c(rep("Black", length(seq(1,starcord1[1]-1)) ), rep(color_arm2,length(seq(starcord1[1], starcord1[2]))),rep(color_arm1,length(seq(maturecord1[1], maturecord1[2]))), rep("Black",length(seq(maturecord1[2], nchar(folded[[1]][1])-1))) )
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
          firstMIRoverlap<-"None"


          ###### Evaluate DROSHA cleavage!

          if(foldingtable_2[firstmirnanuc5p,"dots"]!="."){#if 1st one, has a complementary
            Compl_to_firstmirnanuc5p <- Findmatchingupstream(firstmirnanuc5p) ## Get complementary to first

            if(!is.null(Compl_to_firstmirnanuc5p)){
              if(foldingtable_2[Compl_to_firstmirnanuc5p,"color"]!="Black" & foldingtable_2[Compl_to_firstmirnanuc5p,"color"]!=foldingtable_2[firstmirnanuc5p,"color"] ){# If complementary is not same color not black
                firstMIRoverlap<-"In"
                if(foldingtable_2[Compl_to_firstmirnanuc5p+1,"color"]!="Black" & foldingtable_2[Compl_to_firstmirnanuc5p+2,"color"]!="Black" & foldingtable_2[Compl_to_firstmirnanuc5p+3,"color"]=="Black" ){# if only 2 upstream are also colored
                  if(sum(str_count(foldingtable_2[(firstmirnanuc5p-2):(firstmirnanuc5p+1),"dots"], "\\(" ))==4 ){ # if first, second and 2 previous all matched is perfect
                    print("Perfect Drosha")
                    overhang_animal<-rbind(overhang_animal, c("Perfect pri-miRNA (Drosha) cleavage", perfectmatch_animal) )
                    overhang_plant<-rbind(overhang_plant, c("Perfect pri-miRNA cleavage", perfectmatch_plant))
                  }else if(foldingtable_2[Compl_to_firstmirnanuc5p+3,"color"]=="Black" & sum(str_count(foldingtable_2[(firstmirnanuc5p-2):(firstmirnanuc5p+3),"dots"], "\\(" ))>=3 ){#they will at least have 3 complementary (the first+2) including 3rd before, bc somtimes mini bouble
                    print("Acceptable Drosha")
                    overhang_animal<-rbind(overhang_animal, c("Acceptable pri-miRNA (Drosha) cleavage", acceptable_animal) )
                    overhang_plant<-rbind(overhang_plant, c("Acceptable pri-miRNA cleavage", acceptable_plant))
                  }else if ( sum(str_count(foldingtable_2[(firstmirnanuc5p-2):(firstmirnanuc5p+3),"dots"], "\\(" ))>=3 ) {#complementarity less than 3
                    print("Weak Drosha 1")
                    overhang_animal<-rbind(overhang_animal, c("Suspicious pri-miRNA (Drosha) cleavage", suspicious_animal) )
                    overhang_plant<-rbind(overhang_plant, c("Suspicious pri-miRNA cleavage", suspicious_plant))
                  }else{
                    print("Bad Drosha 4")
                    overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                    overhang_plant<- rbind(overhang_plant, c("Bad pri-miRNA cleavage" , 0))
                  }
                }else if( (foldingtable_2[(Compl_to_firstmirnanuc5p+3),"color"]=="Black" |foldingtable_2[(Compl_to_firstmirnanuc5p+4),"color"]=="Black") & foldingtable_2[Compl_to_firstmirnanuc5p+1,"color"]!="Black") {#If only 1 upstream colored, and 1st is paired
                  print("Weak Drosha 2")
                  overhang_animal<-rbind(overhang_animal, c("Suspicious pri-miRNA (Drosha) cleavage", suspicious_animal) )
                  overhang_plant<-rbind(overhang_plant, c("Suspicious pri-miRNA cleavage", suspicious_plant))
                }else{
                  print("Bad Drosha 1")
                  overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                  overhang_plant<-rbind(overhang_plant, c("Bad pri-miRNA cleavage" , 0))
                }# none upstream colored
                ### if the 1st nuc, complement balck   firstMIRoverlap<-"None"
              }else if(foldingtable_2[Compl_to_firstmirnanuc5p,"color"]=="Black" ) { #if complement of the first is a black
                firstMIRoverlap<-"Out"
                print("Bad Drosha 3p must hang")
                overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha) cleavage, 3p must hang", 0) )
                overhang_plant<- rbind(overhang_plant, c("Bad pri-miRNA cleavage, 3p must hang" , 0))
              }else{ #complement is same color super bad
                print("Bad Drosha 2")
                overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                overhang_plant<- rbind(overhang_plant, c("Bad pri-miRNA cleavage" , 0))
              }
            }else{
              overhang_animal<-rbind(overhang_animal, c("Uncommon structure", 0) )
              overhang_plant<-rbind(overhang_plant, c("Uncommon cleavage" , 0))
            }

          }else if( !is.null(Findmatchingupstream(firstmirnanuc5p+1))){## if first no, pair, but second yes
            if(foldingtable_2[firstmirnanuc5p+1,"dots"]!="." & foldingtable_2[Findmatchingupstream(firstmirnanuc5p+1),"color"]!="Black" & foldingtable_2[Findmatchingupstream(firstmirnanuc5p+1),"color"]!=foldingtable_2[firstmirnanuc5p,"color"]  ){# If  first doesnt have coplementary, check 2nd
              firstMIRoverlap<-"In"
              print("Could still be good")
              overhang_animal<-rbind(overhang_animal, c("Acceptable pri-miRNA (Drosha) cleavage", acceptable_animal) )
              overhang_plant<-rbind(overhang_plant, c("CAcceptable pri-miRNA cleavaeg" , acceptable_animal))
            }else{ #first 2 ones no complement
              print("Bad Drosha 3")
              overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha)cleavage", 0) )
              overhang_plant<-rbind(overhang_plant, c("Bad pri-miRNA cleavage" , 0))
            }
          }else{ #first 2 ones no complement
            print("Bad Drosha 4")
            overhang_animal<-rbind(overhang_animal, c("Bad pri-miRNA (Drosha) cleavage", 0) )
            overhang_plant<-rbind(overhang_plant, c("Bad pri-miRNA cleavage" , 0))
          }

          print(firstMIRoverlap)

          ###### Evaluate DICER cleavage!
          secondMIRoverlap<-"None"
          if (overlap>=-2){# if mature and star overlap, is bad!!
            print("Overlap mature and star")
            overhang2_animal<-rbind(overhang2_animal,c("Bad, No Loop", Nolooppenalty))
            overhang2_plant<-rbind(overhang2_plant,c("Bad, No Loop", Nolooppenalty))
          }else{
            firstcolor<-foldingtable[min(which(foldingtable$color!="Black")),"color"]
            firstmirnaLastnuc=max(which(foldingtable$color==firstcolor))


            if(foldingtable_2[firstmirnaLastnuc,"dots"]!="."){#if 1st one, has a complementary
              Compl_to_firstmirnaLastnuc<-Findmatchingupstream(firstmirnaLastnuc)

              if(!is.null(Compl_to_firstmirnaLastnuc)){

                if(foldingtable_2[Compl_to_firstmirnaLastnuc,"color"]=="Black" ) { #if complement of the first is a black
                  secondMIRoverlap<-"Out"
                  if( foldingtable_2[(Compl_to_firstmirnaLastnuc+1),"color"]=="Black" & foldingtable_2[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black"  & sum(str_count(foldingtable_2[(Compl_to_firstmirnaLastnuc):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))==4  ){### Last one is complement to black (and it is not hanging same arm both times)
                    # if hanging 2 and good matches
                    print("Perfect Dicer 2")
                    overhang2_animal<-rbind(overhang2_animal, c("Perfect loop (Dicer) cleavage", perfectmatch_animal) )
                    overhang2_plant<-rbind(overhang2_plant, c("Perfect loop cleavage", perfectmatch_plant))
                  }else if(  foldingtable_2[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black" & sum(str_count(foldingtable_2[(Compl_to_firstmirnaLastnuc-1):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))>=3){
                    print("Acceptable Dicer 1")
                    overhang2_animal<-rbind(overhang2_animal, c("Acceptable loop (Dicer) cleavage", acceptable_animal) )
                    overhang2_plant<-rbind(overhang2_plant, c("Acceptable loop cleavage", acceptable_plant))
                  }else if( ((foldingtable_2[(Compl_to_firstmirnaLastnuc+3),"color"]!="Black" | foldingtable_2[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black")) & sum(str_count(foldingtable_2[(Compl_to_firstmirnaLastnuc-1):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))>=3) {#
                    #if hanging 3 and few mismatches
                    print("WEAK")
                    overhang2_animal<-rbind(overhang2_animal, c("Suspicious loop (Dicer) cleavage", suspicious_animal) )
                    overhang2_plant<-rbind(overhang2_plant, c("Suspicious loop cleavage", suspicious_plant))
                  }else{
                    print("Bad Dicer 4")
                    overhang2_animal<-rbind(overhang2_animal, c("Bad loop (Dicer) cleavage", 0) )
                    overhang2_plant<-rbind(overhang2_plant, c("Bad loop cleavage" , 0))
                  }
                }else{
                  print("Bad Dicer 5")
                  overhang2_animal<-rbind(overhang2_animal, c("Bad loop cleavage (Dicer), 3p not hanging", 0) )
                  overhang2_plant<-rbind(overhang2_plant, c("Bad Loop cleavage, 3p not hanging" , 0))
                }
              }else{
                overhang2_animal<-rbind(overhang2_animal, c("Uncommon structure", 0) )
                overhang2_plant<-rbind(overhang2_plant, c("Uncommon cleavage" , 0))}

            }else if(foldingtable_2[(firstmirnaLastnuc-1),"dots"]!="."){# if 2nd has a pair
              Compl_to_firstmirnaPENULTtnuc<-Findmatchingupstream(firstmirnaLastnuc-1)
              if (is.null(Findmatchingupstream(firstmirnaLastnuc-1))){
                print("Not normal structure")
                overhang2_animal<-rbind(overhang2_animal, c("Uncommon structure", 0) )
                overhang2_plant<-rbind(overhang2_plant, c("Uncommon cleavage" , 0))
              }else if( foldingtable_2[(Compl_to_firstmirnaPENULTtnuc),"color"]=="Black" &   foldingtable_2[(Compl_to_firstmirnaPENULTtnuc+1),"color"]!="Black" ){
                print("Acceptable Dicer Loop 1b")
                overhang2_animal<-rbind(overhang2_animal, c("Acceptable loop (Dicer) cleavage", acceptable_animal) )
                overhang2_plant<-rbind(overhang2_plant, c("Acceptable loop cleavage", acceptable_plant))

              }else{
                print("Bad Dicer 1b")
                overhang2_animal<-rbind(overhang2_animal, c("Bad loop (Dicer) cleavage", 0) )
                overhang2_plant<-rbind(overhang2_plant, c("Bad loop cleavage" , 0))
              }

            }else{ #complement is same color super bad
              print("Bad Dicer 3")
              overhang2_animal<-rbind(overhang2_animal, c("Bad loop (Dicer) cleavage", 0) )
              overhang2_plant<-rbind(overhang2_plant, c("Bad loop cleavage" , 0))}

          }
        }else{
          print("To many Ns in  precursor!")
          overhang_animal<-rbind(overhang_animal, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang_plant<-rbind(overhang_plant, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang2_animal<-rbind(overhang2_animal, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang2_plant<-rbind(overhang2_plant, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))

        }# if more than 20 Ns in precursor




        ########################check adjust structure###############################
        if (!str_count(mirnadf$precursor[i], "N") > 20){## if miRNA precursorcontains > 20 Nt
          ##if 5P mature
          if(maturecord1_adjust[1]<starcord1_adjust[1]){# if 5p
            print(paste("mature is 5'",i))

            overlap=maturecord1_adjust[2]-starcord1_adjust[1]+1
            if(overlap<0) {#if mature and star don't overlap (as it should)
              colorvector<-c(rep("Black", length(seq(1,maturecord1_adjust[1]-1)) ), rep(color_arm1,length(seq(maturecord1_adjust[1], maturecord1_adjust[2]))),rep("Black",length(seq(maturecord1_adjust[2]+1, starcord1_adjust[1]-1))),rep(color_arm2,length(seq(starcord1_adjust[1], starcord1_adjust[2]))), rep("Black",length(seq(starcord1_adjust[2], nchar(folded[[1]][1])-1))) )
              foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

              foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
              foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
            }else{# If mature and star overlap (they shouldn't!) don't crash do :
              overlaFLAG=TRUE
              if (overlap>0){# if there is overlap
                colorvector<-c(rep("Black", length(seq(1,maturecord1_adjust[1]-1)) ), rep(color_arm2,length(seq(maturecord1_adjust[1], maturecord1_adjust[2]-overlap))),rep("Orange",length(seq(starcord1_adjust[1], maturecord1_adjust[2]))),rep(color_arm1,length(seq(starcord1_adjust[1]+overlap, starcord1_adjust[2]))), rep("Black",length(seq(starcord1_adjust[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
              if (overlap==0){# if they dont overlap, but there is no loop

                colorvector<-c(rep("Black", length(seq(1,maturecord1_adjust[1]-1)) ), rep(color_arm1,length(seq(maturecord1_adjust[1], maturecord1_adjust[2]-overlap))),rep(color_arm2,length(seq(starcord1_adjust[1]+overlap, starcord1_adjust[2]))), rep("Black",length(seq(starcord1_adjust[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
            }



          }else{   ##################### Mature 3'
            print(paste("mature is 3'", i))

            overlap=starcord1_adjust[2]+1-maturecord1_adjust[1]

            if(overlap<0){ # if no overlap

              colorvector<-c(rep("Black", length(seq(1,starcord1_adjust[1]-1)) ), rep(color_arm2,length(seq(starcord1_adjust[1], starcord1_adjust[2]))),rep("Black",length(seq(starcord1_adjust[2]+1, maturecord1_adjust[1]-1))),rep(color_arm1,length(seq(maturecord1_adjust[1], maturecord1_adjust[2]))), rep("Black",length(seq(maturecord1_adjust[2], nchar(folded[[1]][1])-1))) )
              foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

              foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
              foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)

            }else{# If mature and star overlap (they shouldn't!) don't crash do :

              if (overlap>0){
                colorvector<-c(rep("Black", length(seq(1,starcord1_adjust[1]-1)) ), rep(color_arm2,length(seq(starcord1_adjust[1], starcord1_adjust[2]-overlap))),rep("Orange",length(seq(starcord1_adjust[1], maturecord1_adjust[2]))),rep(color_arm1,length(seq(maturecord1_adjust[1]+overlap, maturecord1_adjust[2]))), rep("Black",length(seq(maturecord1_adjust[2], nchar(folded[[1]][1])-1))) )
                foldingtable<- data.frame("color"=colorvector, "dots"=str_split(folded[[1]][2], "")[[1]], "seq"=str_split(folded[[1]][1], "")[[1]] ,"openprent"=NA ,"closeprent"=NA)

                foldingtable[foldingtable$dots=="(",]$"openprent"<-seq(1, nrow( foldingtable[foldingtable$dots=="(",]) )
                foldingtable[foldingtable$dots==")",]$"closeprent"<-seq(nrow( foldingtable[foldingtable$dots==")",]) , 1)
              }
              if (overlap==0){# if overlap is zero but also no gap
                colorvector<-c(rep("Black", length(seq(1,starcord1_adjust[1]-1)) ), rep(color_arm2,length(seq(starcord1_adjust[1], starcord1_adjust[2]))),rep(color_arm1,length(seq(maturecord1_adjust[1], maturecord1_adjust[2]))), rep("Black",length(seq(maturecord1_adjust[2], nchar(folded[[1]][1])-1))) )
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
          foldingtable_2_adjust<-foldingtable
          foldingtable_2_adjust$dotnumeric<-0
          foldingtable_2_adjust[foldingtable_2_adjust$dots=="(",]$dotnumeric<-(-1)
          foldingtable_2_adjust[foldingtable_2_adjust$dots==")",]$dotnumeric<-1

          foldingtable_2_adjust[firstmirnanuc5p,]
          foldingtable_2_adjust

          #Findmatchingupstream(firstmirnanuc5p+1)

          ## find the matching base!
          Findmatchingupstream<-function(querynt){
            for(pos in 1:(nrow(foldingtable_2_adjust)-querynt)){
              if(sum(foldingtable_2_adjust$dotnumeric[(querynt+1) : (querynt+pos) ]) == 1){
                complement=pos+querynt
                return(complement)
                break()
              }
            }
            return(NULL)
          }

          ## I have the 1st nt of the miRNA located at the 5'
          firstMIRoverlap<-"None"


          ###### Evaluate DROSHA cleavage

          if(foldingtable_2_adjust[firstmirnanuc5p,"dots"]!="."){#if 1st one, has a complementary
            Compl_to_firstmirnanuc5p <- Findmatchingupstream(firstmirnanuc5p) ## Get complementary to first

            if(!is.null(Compl_to_firstmirnanuc5p)){
              if(foldingtable_2_adjust[Compl_to_firstmirnanuc5p,"color"]!="Black" & foldingtable_2_adjust[Compl_to_firstmirnanuc5p,"color"]!=foldingtable_2_adjust[firstmirnanuc5p,"color"] ){# If complementary is not same color not black
                firstMIRoverlap<-"In"
                if(foldingtable_2_adjust[Compl_to_firstmirnanuc5p+1,"color"]!="Black" & foldingtable_2_adjust[Compl_to_firstmirnanuc5p+2,"color"]!="Black" & foldingtable_2_adjust[Compl_to_firstmirnanuc5p+3,"color"]=="Black" ){# if only 2 upstream are also colored
                  if(sum(str_count(foldingtable_2_adjust[(firstmirnanuc5p-2):(firstmirnanuc5p+1),"dots"], "\\(" ))==4 ){ # if first, second and 2 previous all matched is perfect
                    print("Perfect Drosha")
                    overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Perfect pri-miRNA (Drosha) cleavage", perfectmatch_animal) )
                    overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Perfect pri-miRNA cleavage", perfectmatch_plant))
                  }else if(foldingtable_2_adjust[Compl_to_firstmirnanuc5p+3,"color"]=="Black" & sum(str_count(foldingtable_2_adjust[(firstmirnanuc5p-2):(firstmirnanuc5p+3),"dots"], "\\(" ))>=3 ){#they will at least have 3 complementary (the first+2) including 3rd before, bc somtimes mini bouble
                    print("Acceptable Drosha")
                    overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Acceptable pri-miRNA (Drosha) cleavage", acceptable_animal) )
                    overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Acceptable pri-miRNA cleavage", acceptable_plant))
                  }else if ( sum(str_count(foldingtable_2_adjust[(firstmirnanuc5p-2):(firstmirnanuc5p+3),"dots"], "\\(" ))>=3 ) {#complementarity less than 3
                    print("Weak Drosha 1")
                    overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Suspicious pri-miRNA (Drosha) cleavage", suspicious_animal) )
                    overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Suspicious pri-miRNA cleavage", suspicious_plant))
                  }else{
                    print("Bad Drosha 4")
                    overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                    overhang_plant_adjust<- rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage" , 0))
                  }
                }else if( (foldingtable_2_adjust[(Compl_to_firstmirnanuc5p+3),"color"]=="Black" |foldingtable_2_adjust[(Compl_to_firstmirnanuc5p+4),"color"]=="Black") & foldingtable_2_adjust[Compl_to_firstmirnanuc5p+1,"color"]!="Black") {#If only 1 upstream colored, and 1st is paired
                  print("Weak Drosha 2")
                  overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Suspicious pri-miRNA (Drosha) cleavage", suspicious_animal) )
                  overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Suspicious pri-miRNA cleavage", suspicious_plant))
                }else{
                  print("Bad Drosha 1")
                  overhang_animal_adjust <- rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                  overhang_plant_adjust <- rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage" , 0))
                }# none upstream colored
                ### if the 1st nuc, complement balck   firstMIRoverlap<-"None"
              }else if(foldingtable_2_adjust[Compl_to_firstmirnanuc5p,"color"]=="Black" ) { #if complement of the first is a black
                firstMIRoverlap<-"Out"
                print("Bad Drosha 3p must hang")
                overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha) cleavage, 3p must hang", 0) )
                overhang_plant_adjust<- rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage, 3p must hang" , 0))
              }else{ #complement is same color super bad
                print("Bad Drosha 2")
                overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha) cleavage", 0) )
                overhang_plant_adjust<- rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage" , 0))
              }
            }else{
              overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Uncommon structure", 0) )
              overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Uncommon cleavage" , 0))
            }

          }else if( !is.null(Findmatchingupstream(firstmirnanuc5p+1))){## if first no, pair, but second yes
            if(foldingtable_2_adjust[firstmirnanuc5p+1,"dots"]!="." & foldingtable_2_adjust[Findmatchingupstream(firstmirnanuc5p+1),"color"]!="Black" & foldingtable_2_adjust[Findmatchingupstream(firstmirnanuc5p+1),"color"]!=foldingtable_2_adjust[firstmirnanuc5p,"color"]  ){# If  first doesnt have coplementary, check 2nd
              firstMIRoverlap<-"In"
              print("Could still be good")
              overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Acceptable pri-miRNA (Drosha) cleavage", acceptable_animal) )
              overhang_plant_adjust<-rbind(overhang_plant_adjust, c("CAcceptable pri-miRNA cleavaeg" , acceptable_plant))
            }else{ #first 2 ones no complement
              print("Bad Drosha 3")
              overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha)cleavage", 0) )
              overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage" , 0))
            }
          }else{ #first 2 ones no complement
            print("Bad Drosha 4")
            overhang_animal_adjust<-rbind(overhang_animal_adjust, c("Bad pri-miRNA (Drosha) cleavage", 0) )
            overhang_plant_adjust<-rbind(overhang_plant_adjust, c("Bad pri-miRNA cleavage" , 0))
          }

          print(firstMIRoverlap)

          ###### Evaluate DICER cleavage!
          secondMIRoverlap<-"None"
          if (overlap>=-2){# if mature and star overlap, is bad!!
            print("Overlap mature and star")
            overhang2_animal_adjust<-rbind(overhang2_animal_adjust,c("Bad, No Loop", Nolooppenalty))
            overhang2_plant_adjust<-rbind(overhang2_plant_adjust,c("Bad, No Loop", Nolooppenalty))
          }else{
            firstcolor<-foldingtable[min(which(foldingtable$color!="Black")),"color"]
            firstmirnaLastnuc=max(which(foldingtable$color==firstcolor))


            if(foldingtable_2_adjust[firstmirnaLastnuc,"dots"]!="."){#if 1st one, has a complementary
              Compl_to_firstmirnaLastnuc<-Findmatchingupstream(firstmirnaLastnuc)

              if(!is.null(Compl_to_firstmirnaLastnuc)){

                if(foldingtable_2_adjust[Compl_to_firstmirnaLastnuc,"color"]=="Black" ) { #if complement of the first is a black
                  secondMIRoverlap<-"Out"
                  if( foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc+1),"color"]=="Black" & foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black"  & sum(str_count(foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))==4  ){### Last one is complement to black (and it is not hanging same arm both times)
                    # if hanging 2 and good matches
                    print("Perfect Dicer 2")
                    overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Perfect loop (Dicer) cleavage", perfectmatch_animal) )
                    overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Perfect loop cleavage", perfectmatch_plant))
                  }else if(  foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black" & sum(str_count(foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc-1):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))>=3){
                    print("Acceptable Dicer 1")
                    overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Acceptable loop (Dicer) cleavage", acceptable_animal) )
                    overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Acceptable loop cleavage", acceptable_plant))
                  }else if( ((foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc+3),"color"]!="Black" | foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc+2),"color"]!="Black")) & sum(str_count(foldingtable_2_adjust[(Compl_to_firstmirnaLastnuc-1):(Compl_to_firstmirnaLastnuc+3),"dots"], "\\)" ))>=3) {#
                    #if hanging 3 and few mismatches
                    print("WEAK")
                    overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Suspicious loop (Dicer) cleavage", suspicious_animal) )
                    overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Suspicious loop cleavage", suspicious_plant))
                  }else{
                    print("Bad Dicer 4")
                    overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Bad loop (Dicer) cleavage", 0) )
                    overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Bad loop cleavage" , 0))
                  }
                }else{
                  print("Bad Dicer 5")
                  overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Bad loop cleavage (Dicer), 3p not hanging", 0) )
                  overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Bad Loop cleavage, 3p not hanging" , 0))
                }
              }else{
                overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Uncommon structure", 0) )
                overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Uncommon cleavage" , 0))}

            }else if(foldingtable_2_adjust[(firstmirnaLastnuc-1),"dots"]!="."){# if 2nd has a pair
              Compl_to_firstmirnaPENULTtnuc<-Findmatchingupstream(firstmirnaLastnuc-1)
              if (is.null(Findmatchingupstream(firstmirnaLastnuc-1))){
                print("Not normal structure")
                overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Uncommon structure", 0) )
                overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Uncommon cleavage" , 0))
              }else if( foldingtable_2_adjust[(Compl_to_firstmirnaPENULTtnuc),"color"]=="Black" &   foldingtable_2_adjust[(Compl_to_firstmirnaPENULTtnuc+1),"color"]!="Black" ){
                print("Acceptable Dicer Loop 1b")
                overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Acceptable loop (Dicer) cleavage", acceptable_animal) )
                overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Acceptable loop cleavage", acceptable_plant))

              }else{
                print("Bad Dicer 1b")
                overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Bad loop (Dicer) cleavage", 0) )
                overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Bad loop cleavage" , 0))
              }

            }else{ #complement is same color super bad
              print("Bad Dicer 3")
              overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c("Bad loop (Dicer) cleavage", 0) )
              overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c("Bad loop cleavage" , 0))}

          }
        }else{
          print("To many Ns in  precursor!")
          overhang_animal_adjust<-rbind(overhang_animal_adjust, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang_plant<-rbind(overhang_plant, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang2_animal_adjust<-rbind(overhang2_animal_adjust, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))
          overhang2_plant_adjust<-rbind(overhang2_plant_adjust, c(paste("Too many Ns in prec", str_count(mirnadf$precursor[i], "N")), 0))

        }# if more than 20 Ns in precursor


        ##################caculate whether we get at least 16nts complementary
        countComp1 <- 0 # check the user-provided annotation
        countComp2 <- 0 # check the adjust annotation
        for (j in maturecord1[1] : maturecord1[2]) {
          if ( coord$bound[j] %in% starcord1[1] : starcord1[2] ) {
            countComp1 = countComp1 + 1
          }
        }
        for (j in maturecord1_adjust[1] : maturecord1_adjust[2]) {
          if ( coord$bound[j] %in% starcord1_adjust[1] : starcord1_adjust[2] ) {
            countComp2 = countComp2 + 1
          }
        }

        if (countComp1 < 16) {
          overhang2_animal[i,2] <- -5
          overhang2_plant[i,2] <- -5

        }

        if (countComp2 < 16) {
          overhang2_animal_adjust[i, 2] <- -5
          overhang2_plant_adjust[i, 2] <- -5

        }


        ###############creat a final list contain structure information ####################


        if (specie == "Animal" ) {
          if ((as.numeric(overhang2_animal[i,2]) + as.numeric(overhang_animal[i, 2])) < (as.numeric(overhang2_animal_adjust[i, 2])+ as.numeric(overhang_animal_adjust[i, 2]))) {
            print(((as.numeric(overhang2_animal[i,2]) + as.numeric(overhang_animal[i, 2])) < (as.numeric(overhang2_animal_adjust[i, 2])+ as.numeric(overhang_animal_adjust[i, 2]))))
            levels(mirnadf$mature) <- c(levels(mirnadf$mature), adjustMatureSequence [i])
            mirnadf$mature[i] = adjustMatureSequence [i]
            levels(mirnadf$star) <- c(levels(mirnadf$star), adjustStarSequence [i])
            mirnadf$star [i] = adjustStarSequence [i]
            overhang_animal[i,] <- overhang_animal_adjust[i,]
            overhang2_animal [i,]<- overhang2_animal_adjust[i,]
            finalStarPosition[[i]] <- adjustStarPosition[[i]]
            finalMaturePosition[[i]] <- adjustMaturePosition[[i]]
            foldingtable_2 <- foldingtable_2_adjust

          } else {
            finalStarPosition [[i]] <- starcord1
            finalMaturePosition[[i]] <- maturecord1
          }
        } else {

          if ((as.numeric(overhang2_plant[i,2]) + as.numeric(overhang_plant[i, 2])) < (as.numeric(overhang2_plant_adjust[i, 2])+ as.numeric(overhang_plant_adjust[i, 2]))) {
            print(((as.numeric(overhang2_plant[i,2]) + as.numeric(overhang_plant[i, 2])) < (as.numeric(overhang2_plant_adjust[i, 2])+ as.numeric(overhang_plant_adjust[i, 2]))))
            levels(mirnadf$mature) <- c(levels(mirnadf$mature), adjustMatureSequence [i])
            mirnadf$mature[i] = adjustMatureSequence [i]
            levels(mirnadf$star) <- c(levels(mirnadf$star), adjustStarSequence [i])
            mirnadf$star [i] = adjustStarSequence [i]
            overhang_plant[i,] <- overhang_plant_adjust[i,]
            overhang2_plant[i,]<- overhang2_plant_adjust[i,]
            finalStarPosition[[i]] <- adjustStarPosition[[i]]
            finalMaturePosition[[i]] <- adjustMaturePosition[[i]]
            foldingtable_2 <- foldingtable_2_adjust
          } else {
            finalStarPosition [[i]] <- starcord1
            finalMaturePosition[[i]] <- maturecord1
          }
        }


        ## Check the mature/star sequence length


        penaltyscore1 [i]<-ifelse((max(finalStarPosition[[i]]) - min(finalStarPosition[[i]]) + 1) < 18 |  (max(finalStarPosition[[i]]) - min(finalStarPosition[[i]]) + 1) > 26, penalty_length, 0 )
        penaltyscore2 [i]<-ifelse((max(finalMaturePosition[[i]]) - min(finalMaturePosition[[i]]) + 1) < 18 |  (max(finalMaturePosition[[i]]) - min(finalMaturePosition[[i]]) + 1) > 26, penalty_length, 0 )
        penaltyscorelength_animal<<-penaltyscore2+penaltyscore1

        ## later I would add plants



        ###Lets adjust image parameters depending on length
        if(nchar(folded[1,]) > 180){

          jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
          par(mar=c(0.1,0.1,0,0.1))
          RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
                  seqcols=c(color_arm1,color_arm2),labTF=FALSE,
                  pointSize = 1, lineWd = 1, nt=T,
                  dp=1, tsize=0.5)
          dev.off()
        }else if( nchar(folded[1,]) > 100){
          jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
          par(mar=c(0.1,0.1,0,0.1))
          RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
                  seqcols=c(color_arm1,color_arm2),labTF=FALSE,
                  pointSize = 2, lineWd = 1, nt=T,
                  dp=1, tsize=0.8)
          dev.off()
        }else{
          jpeg(filename = paste( "www/images/", mirnadf$ID[i],"_fold.jpg",sep=''),quality=100, width = 2000, height = 2000, units = "px",res =300 )
          par(mar=c(0.1,0.1,0,0.1))
          RNAPlot(coord,hl=c(as.character(RNAString(DNAString(mirnadf$mature[i]))), as.character(RNAString(DNAString(mirnadf$star[i])))),#, main=mirnadf$ID[i]
                  seqcols=c(color_arm1,color_arm2),labTF=FALSE,
                  pointSize = 2.3, lineWd = 1, nt=T,
                  dp=1, tsize=1)
          dev.off()
        }
        folded_globe [[i]] <- folded

        ##################I should use it##############################


        incProgress(1/q, detail = paste("Prec", i, "of", q))


      }# close loop for each prec
      finalMaturePosition <<- finalMaturePosition
      finalStarPosition <<- finalStarPosition
    })# close section folding

    print(paste("Good here!! 1", i))
    foldingFigs<-paste("<img src=\"images/",mirnadf$ID,"_fold.jpg\" width=\"500\ height=\"180\"></img>", sep="")



    ## scores animal / plant#######################################################
    if (specie=="Animal"){
      print("Making dataframe animal")
      print(overhang_animal)
      print(overhang2_animal)
      mirnadf_folding<-cbind(mirnadf,foldingFigs ,"pri-miRNA Cleavage"=overhang_animal[,1],"Loop Cleavage"=overhang2_animal[,1] )
      output$mirnaSeqswithplots <-  DT::renderDataTable({ mirnadf_folding[,c(1,2,3,7,8,9)]},  escape = FALSE )
      overhangs_score_animal<<- as.numeric(overhang_animal[,2])+as.numeric(overhang2_animal[,2] )
    }else{
      print("Making dataframe plant")
      mirnadf_folding<-cbind(mirnadf,foldingFigs ,"pri-miRNA Cleavage"=overhang_plant[,1],"Loop Cleavage"=overhang2_plant[,1] )
      output$mirnaSeqswithplots <-  DT::renderDataTable({ mirnadf_folding[,c(1,2,3,7,8,9)]},  escape = FALSE )
      overhangs_score_plant<<- as.numeric(overhang_plant[,2])+as.numeric(overhang2_plant[,2] )
    }

    foldingFigs<<-foldingFigs

    values$successStep2<-TRUE
    print(  "values$successStep2==TRUE")
    showNotification("Done. Check RNA folding Tab", type= "message")

  })## clsose button fold


  #########################final expression part##################################


  observeEvent(input$ButtonExp, {
    shiny::validate(
      need( values$successStep1==TRUE, message = ('Missing succesful Step 1')),
      errorClass =  showNotification("Missing succesful Step 1", type= "error")
    )
    shiny::validate(
      need( values$successStepAdjust==TRUE, message = ('Missing succesful Step 3: fold seqs ')),
      errorClass =  showNotification("Missing succesful Step 3:fold seqs", type= "error")
    )
    Score_expression_animal<-rep(0, nrow(mirnadf)) # set all scores to zero
    Score_expression_plant<-rep(0, nrow(mirnadf)) # set all scores to zero
    withProgress(message = 'calculating Expression...', value = 0, {

      q <- nrow(mirnadf)
      ### Make it fast for trials

      matureCounts <- c()
      starCounts <- c()
      loopCounts <- c()
      loopReason <- rep(0, nrow(mirnadf))
      flankReason <- rep(0, nrow(mirnadf))
      homology <- rep(0, nrow(mirnadf))
      reason <- c()

      for(i in 1 : nrow(mirnadf) ) {
        incProgress(1/q, detail = paste("Plot", i, "of", q))
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
        plot<-barplot(toplot, axes=TRUE, ylab="Number of Reads", main=mirname, col=c("blue4", "grey"), border=c("blue4","grey"),beside=T)
        # mtext(at = plot, text = seq,col="black", side = 1,  line = 0, cex=1)

        stardf_extrabases<-GRanges(stardf)
        matdf_extrabases<-GRanges(matdf)



        star_i <- finalStarPosition[[i]][1]
        star_e <- finalStarPosition[[i]][2]
        matureseq_i <- finalMaturePosition[[i]][1]
        matureseq_e <- finalMaturePosition[[i]][2]


        mtext(at = plot[1,c(1:matureseq_i,matureseq_e:length(seq))], text =seq[c(1:matureseq_i,matureseq_e:length(seq))] ,col="black", side = 1,  line = 0, cex=1)
        mtext(at = plot[1,star_i:star_e], text = seq[star_i:star_e],col=color_arm2, side = 1,  line = 0, cex=1, font=( face=2))
        mtext(at = plot[1,matureseq_i:matureseq_e], text = seq[matureseq_i:matureseq_e],col=color_arm1, side = 1,  line = 0, cex=1, font=( face=2))
        dev.off()


        ##############check the expression###################
        starReads <- mean(selectedRange_coverage[star_i:star_e])
        matureReads <- mean(selectedRange_coverage[matureseq_i:matureseq_e])

        if (star_i >= matureseq_e) {
          loopReads <- mean(selectedRange_coverage[(matureseq_e+1 + 1):(star_i-3)])
          flankReadsLeft <- mean(selectedRange_coverage[(star_e+1):length(selectedRange_coverage)])
          flankReadsright <- mean(selectedRange_coverage[1:(matureseq_i-1)])
        } else {
          loopReads <- mean(selectedRange_coverage[(star_e+ 1 + 1):(matureseq_i-3)])
          flankReadsLeft <- mean(selectedRange_coverage[(matureseq_e+1):length(selectedRange_coverage)])
          flankReadsright <- mean(selectedRange_coverage[1:(star_i-1)])
        }

        if (starReads > matureReads) {
          tempReads <-  matureReads
          matureReads <- starReads
          starReads <- tempReads

          checkMatureReadsRight <- mean(selectedRange_coverage[(matureseq_e-1):matureseq_e])
          checkMatureReadsLeft <- mean(selectedRange_coverage[matureseq_i:(matureseq_i+1)])
          checkStarReadsRight <- mean(selectedRange_coverage[(star_e-1):star_e])
          checkStarReadsLeft <- mean(selectedRange_coverage[star_i:(star_i+1)])

          tempReads <- checkMatureReadsRight
          checkMatureReadsRight <- checkStarReadsRight
          checkStarReadsRight <- tempReads

          tempReads <- checkMatureReadsLeft
          checkMatureReadsLeft <- checkStarReadsLeft
          checkStarReadsLeft <- tempReads

        } else {

          checkMatureReadsRight <- mean(selectedRange_coverage[(matureseq_e-1):matureseq_e])
          checkMatureReadsLeft <- mean(selectedRange_coverage[matureseq_i:(matureseq_i+1)])
          checkStarReadsRight <- mean(selectedRange_coverage[(star_e-1):star_e])
          checkStarReadsLeft <- mean(selectedRange_coverage[star_i:(star_i+1)])

        }

        ## make sure mature reads is the largest!!!
        matureCounts[i] <- round(matureReads)
        starCounts[i] <- round(starReads)
        loopCounts[i] <- round(loopReads,2)
        ######check the mountain-like structure#########

        #########Give the score For expression############################
        if (starReads == 0 | matureReads == 0) { # IF no expression evidence, give a punishment!
          scoreExpression <- c (-5)
        } else {
          ## Check the expression and loop reads
          if(starReads > ExpLevel2 & matureReads > ExpLevel2 ) {
            scoreExpression <- c(2.5)
          } else if (starReads > ExpLevel1 & matureReads > ExpLevel1) {
            scoreExpression <- c (2)
          } else if (starReads > ExpLevel1 & matureReads > ExpLevel2) {
            scoreExpression <- c (2)
          } else {
            scoreExpression <- c(0)
          }

          ## check loop Reads
          if ((loopReads / (starReads + matureReads + loopReads) > 0.2 )) {
            loopReason[i] = 1
            scoreExpression = scoreExpression - 2
            print("Too many reads in the loop")
          } else if (loopReads / (starReads + matureReads + loopReads) > 0.1) {
            loopReason[i] = 1
            scoreExpression = scoreExpression - 1.5
            print("Too many reads in the loop")
          } else if (loopReads / (starReads + matureReads + loopReads) > 0.06) {
            loopReason[i] = 1
            scoreExpression = scoreExpression - 1.25
            print("Too many reads in the loop")
          } else if (loopReads / (starReads + matureReads + loopReads) > 0.03) {
            loopReason[i] = 1
            scoreExpression = scoreExpression - 1
            print("Too many reads in the loop")
          } else if (loopReads / (starReads + matureReads + loopReads) > 0.012) {
            print("Too many reads in the loop")
            scoreExpression = scoreExpression - 0.75

          }




          ## Check mountain-like structure

            if ((checkMatureReadsLeft / matureReads < 0.9) | (checkStarReadsLeft / starReads < 0.9)) {
              homology[i] = 1
              scoreExpression = scoreExpression - 0.5
              print("Penalty mountain-like expression 0.5")
            }

          ###check franking expression
          if (flankReadsLeft > flankReadsright) {
            if((flankReadsLeft / (matureReads+starReads) > 0.4)) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 2
              print("Penalty franking expression 2")

            } else if ((flankReadsLeft/ (matureReads+starReads) > 0.2) ){
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 1.5
              print("Penalty franking expression 1.5")

            } else if ((flankReadsLeft/ (matureReads+starReads) > 0.12) ) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 0.75
              print("Penalty franking expression 0.75")

            } else if ((flankReadsLeft/ (matureReads+starReads) > 0.04) ) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 0.25
              print("Penalty franking expression 0.25")
            }
          } else {
            if( (flankReadsright / (matureReads + starReads) > 0.4 )) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 2
              print("Penalty franking expression 2")

            } else if ( (flankReadsright / (matureReads + starReads) > 0.2)){
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 1.5
              print("Penalty franking expression 1.5")

            } else if ((flankReadsright/ (matureReads+starReads) > 0.12 )) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 0.75
              print("Penalty franking expression 0.75")

            } else if ((flankReadsright/ (matureReads+starReads) > 0.04 )) {
              flankReason[i] = 1
              scoreExpression <- scoreExpression - 0.25
              print("Penalty franking expression 0.25")
            }

            if ( loopReads / starReads > 0.5) {
              scoreExpression <- scoreExpression - 1
              print("not enough reads in the star")
            }
          }
        }
        Score_expression_animal[i] <- scoreExpression
        Score_expression_plant[i] <- scoreExpression
      }#close for
      Score_expression_animal <<- Score_expression_animal
      loopCounts <<- loopCounts
      starCounts <<- starCounts
      matureCounts <<- matureCounts
      ### Find whis is mature/star
      counts_Table<<-data.frame(loopCounts, matureCounts,starCounts)
      colnames(counts_Table)<- c("Loop","Mature","Star")

      #################### if input data was 5p / 3p, check whch is the mature/star###
      maturesvector <<- NULL
      maturesvector$Id<-as.character(mirnadf$ID)
      if (input$matureorarm == "arm5p3p"){
        maturesvector$whoismature<-ifelse(counts_Table$Mature>counts_Table$Star, "5p is Mature", "3p is Mature" )
        maturesvector$matureis<-ifelse(counts_Table$Mature>counts_Table$Star, "5p", "3p" )
      }
      maturesvector<<-as.data.frame(maturesvector)

      #################################################################################
      ExpressionPlot= paste("<img src=\"plots/",mirnadf$ID,".png\" width=\"1000\" height=\"600\"></img>", sep="")
      if (input$matureorarm == "maturestar"){ # if input data was mature/star
        mirnadf_plots<-data.frame("ID"=mirnadf$ID,"Mature seq"=mirnadf$mature,"Star seq"=mirnadf$star, "Reads loop"=counts_Table$Loop, "Reads Mature"=counts_Table$Mature, "Reads Star"=counts_Table$Star)
        output$PLOTS <-  DT::renderDataTable({ mirnadf_plots},  escape = FALSE, selection = 'single' )
        ExpressionPlot<<-ExpressionPlot
      }else{# if was 5p 3p
        mirnadf_plots<-data.frame("ID"=mirnadf$ID,"5P arm seq"=mirnadf$mature,"3P arm seq"=mirnadf$star, "Reads Mature"= maturesvector$matureis, "Reads loop"=counts_Table$Loop, "Reads 5P arm"=counts_Table$Mature, "Reads 3P arm "=counts_Table$Star)
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


    })

  })
















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



        toalignwhole<<-mirNAStoalign[grep(matureasRNA, as.character(mirNAStoalign)), ] ### select identical mirnas

        toalignseed<<-mirNAStoalign[grep(paste("^.",seed1,sep=""), as.character(mirNAStoalign)), ] ### select identical mirnas

        #names(matureasRNA)<-paste0("<div> <span style=\"color:red\"><strong>CANDIDATE </strong> </span></div>")
        names(matureasRNA)<-"CANDIDATE"

        toalign<-RNAStringSet(matureasRNA)

        numberofalignments<-20

        if (length(toalignwhole)>1 ){ # if some are >1 identical hits
          if (length(toalignwhole)<numberofalignments &  (length(toalignseed)+length(toalignwhole))>numberofalignments  ){ # if identical are less than numberofalignments but plus same seed more than numberofalignments
            toalign<-c(toalign,toalignwhole )####
            toalign<-c(toalign, sample(toalignseed,numberofalignments-length(toalignwhole)) )#### we ad same seed until numberofalignments

            if(length(toalignwhole)>5){## if there are between 6 and 20 identicals  (+ >20 similar)
              score2_animal<<-c(score2_animal,Score_conservation_6_20id20_Animal)
              score2_plant<<-c(score2_plant,Score_conservation_6_20id20_Plant)
              conservationtype<<-c(conservationtype,"strong 6")
            }else{
              score2_animal<<-c(score2_animal,Score_conservation_2_5id_Animal)###### if there are between 2-5 identicals (+ >20 similar)
              score2_plant<<-c(score2_plant,Score_conservation_2_5id_Plant)
              conservationtype<<-c(conservationtype,"strong 4")

            }
          }

          if (length(toalignwhole)<numberofalignments & length(toalignwhole)>0 & length(toalignseed)+length(toalignwhole)<numberofalignments  ){
            toalign<-c(toalign, toalignwhole)####if identical + same seed less numberofalignments, we use all
            toalign<-c(toalign, toalignseed)


            if(length(toalignwhole)>5){## if there are between 6 and 20 identicals  (+ <20 similar)
              score2_animal<<-c(score2_animal, Score_conservation_6_20id_Animal)
              score2_plant<<-c(score2_plant, Score_conservation_6_20id_Plant)
              conservationtype<<-c(conservationtype,"strong 5.8")

            }else{
              score2_animal<<-c(score2_animal,Score_conservation_2_5id_Animal)###### if there are between 2-4 identicals (+ <20 similar)
              score2_plant<<-c(score2_plant,Score_conservation_2_5id_Plant)
              conservationtype<<-c(conservationtype,"medium 4.8")
            }

          }
          if(length(toalignwhole)>numberofalignments){#### if more than numberofalignments identicals, we select ALL identical
            toalign<-c(toalign,toalignwhole)
            score2_animal<<-c(score2_animal,Score_conservation_20id_Animal)######
            score2_plant<<-c(score2_plant,Score_conservation_20id_Plant)######
            conservationtype<<-c(conservationtype,"very strong 6.5")
          }

        }else if(length(toalignseed)>1){#### if no identicals but some with same seed
          if (length(toalignseed)<numberofalignments  ){#if less than numberofalignments we use all
            toalign<-c(toalign,toalignseed)
            score2_animal<<-c(score2_animal,Score_conservation_0id_Animal)######
            score2_plant<<-c(score2_plant, Score_conservation_0id_Plant)######
            conservationtype<<-c(conservationtype,"low 2")


          }
          if (length(toalignseed)>numberofalignments  ){#if more same seed than numberofalignments(20), we select ALL
            toalign<- c(toalign,toalignseed)###
            score2_animal<<-c(score2_animal,Score_conservation_0id20_Animal)######
            score2_plant<<-c(score2_plant,Score_conservation_0id20_Plant)######
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

    CalculateScore<<-function(expression, overhang,homology,penaltylen){
      finalscore<<-expression+overhang+homology+penaltylen

      return(finalscore)
    }

    if(specie=="Animal"){
      Score<-Score_expression_animal+overhangs_score_animal+Score_homology_animal+penaltyscorelength_animal
      mirnadf_integrated_globalcopie<<-cbind(mirnadf,counts_Table,conservationtype,Score, Score_expression_animal, overhangs_score_animal, Score_homology_animal, penaltyscorelength_animal)

    }else{
      Score<-Score_expression_plant+overhangs_score_plant+Score_homology_plant
      mirnadf_integrated_globalcopie<<-cbind(mirnadf,conservationtype,Score, Score_expression_plant, overhangs_score_plant, Score_homology_plant)

    }

    mirnadf_integrated<-cbind(mirnadf,counts_Table,conservationtype,Score)


    print("Integrating")
    #### Threshold
    scorethreshold <<- input$threshold#3.75

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
    ## Download PDF report
    if (specie == "Animal") {
      trueMiRNA <- toreturnall[toreturnall$Score > scorethreshold | (toreturnall$Score_expression_animal > 2 & toreturnall$overhangs_score_animal > 1.25),]
    } else {
      trueMiRNA <- toreturnall[toreturnall$Score > scorethreshold | (toreturnall$Score_expression_plant > 2 & toreturnall$overhangs_score_plant > 1.25),]
    }

    observeEvent (input$report, {
      withProgress(message = 'Creating report...', value = 0, {
        n <- length(as.character(trueMiRNA$ID))
        for (i in as.character(trueMiRNA$ID)) {
          render("./report.Rmd", output_file = paste0 ("./report/report_", i, ".pdf"),
                 params = list(new_title = paste ("report of ", i)))
          incProgress(1/n, detail = paste("report for", i))

        }
      })
      showNotification("Done. Your report is in the \"./mirPlot/report\"file", type= "message")

    })

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
        stat_summary(fun.y=mean, colour=color_arm1, geom="point") +
        geom_boxplot (aes(fill=Score), alpha=.5, width=1, position = position_dodge(width = 1),  outlier.colour = "dark gray", outlier.size = 1)+
        ggtitle("Scores per class") + theme(legend.position="none",axis.text=element_text(size=13), axis.title=element_text(size=15))})

    output$ScoreBoxplots2<- renderPlot({

      ggplot(Scoresboxplotdata, aes(x = class , y = log(Score))) +
        stat_summary(fun.y=mean, colour=color_arm1, geom="point") +
        geom_boxplot (aes(fill=Score), alpha=.5, width=1, position = position_dodge(width = 1),  outlier.colour = "dark gray", outlier.size = 1)+
        ggtitle("Log(Scores per class)") + theme(legend.position="none",axis.text=element_text(size=13), axis.title=element_text(size=15))
    })




  })# close button integration

}#close server







