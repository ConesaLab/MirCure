# Define UI for dataset viewer app ----
library("shinythemes")
ui <- navbarPage("MirCure",theme = shinytheme("flatly"),

                 tabPanel(title='Load Data',
                          sidebarLayout(
                            sidebarPanel(

                              radioButtons("specie", "What kind of organism?",
                                           c("Animal" = "Animal","Plant" = "Plant" )),



                              radioButtons("Genome", "Genome location:",
                                           c("Select available genome" = "Select","Upload a genome" = "Upload" )),

                              conditionalPanel(
                                condition = "input.Genome == 'Upload'",
                                fileInput("genome0", label="Upload your genome file",  multiple = FALSE, accept = c(".fa",".fasta"),width="60%")
                              ) ,
                              conditionalPanel(
                                condition = "input.Genome == 'Select'",
                                selectInput("genome1", "Choose Genome from List",choices = c("NULL",list.files('data/genomes/')), selected = NULL ,  multiple = FALSE,width="60%")
                              ) ,

                              tags$hr(),
                              # Input: Select a file ----
                              fileInput("precs", "Precursors (.gff3)",
                                        multiple = FALSE,
                                        accept = c(".gff3"),width="60%"),
                              # Input: Select a file ----
                              fileInput("mature", "Mature or 5' arm (.gff3)",
                                        multiple = FALSE,
                                        accept = c(".gff3"),width="60%"),
                              # Input: Select a file ----
                              fileInput("star", "Star or 3' arm (.gff3)",
                                        multiple = FALSE,
                                        accept = c(".gff3"),width="60%"),
                              # Horizontal line ----
                              #tags$hr(),
                              radioButtons("matureorarm", "Are you uploading mature/star sequences or 5p/3p ?",
                                           c("5p/3p" = "arm5p3p" ,"Mature/Star" = "maturestar"), selected = "arm5p3p"),


                              numericInput("extrabases", "Precursor flanking bases to retrieve:", 11, min = 0, max = 20,width="30%"),

                              actionButton("ButtonSeqs", "1- Get sequences"),
                              # Horizontal line ----

                              radioButtons("Bammenu", "Type:",
                                           c("Select available BAM" = "Selectbam","Upload your BAM" = "Uploadbam" )),

                              conditionalPanel(
                                condition = "input.Bammenu == 'Uploadbam'",
                                fileInput("bam", "BAM file (.BAM)", multiple = FALSE,accept = c(".bam"),width="60%")
                              ) ,
                              conditionalPanel(
                                condition = "input.Bammenu == 'Selectbam'",
                                selectInput("bam1", "Choose BAM from List",choices = c("NULL",list.files('data/bamfiles/')), selected = NULL ,  multiple = FALSE,width="60%")
                              ) ,
                              actionButton("ButtonCheck", "2- Adjust Structure"),

                              tags$hr(),
                              actionButton("ButtonFold", "3- Fold sequences (might take a long time)"),
                              # Horizontal line ----
                              tags$hr(),


                              # Input: Select a file ----

                              #        tags$hr(),
                              actionButton("ButtonExp", "4- Calculate expression"),
                              tags$hr(),
                              actionButton("Homology", "5- Conservation"),
                              tags$hr(),
                              numericInput("threshold", "Score threshold:", 3.75, min = -100, max = 100, step=0.05,width="30%"),
                              actionButton("ButtonIntegrate", "6- Intergrate")
                              #tags$hr(),
                              #bookmarkButton()

                            ),

                            # Main panel for displaying outputs ----
                            mainPanel(
                              h1("MirCure"),

                              fluidRow(
                                column(6,
                                       h4("0 - Load candidate miRNA annotations"),
                                       tags$li(" Upload the genome file (.fa)."),
                                       tags$li(" The miRNA precursor annotations* (.gff3)."),
                                       tags$li(" The miRNA mature annotations* (.gff3)."),
                                       tags$li(" The miRNA star annotations* (.gff3)."),
                                       h6("* All gff3 entries must be in the same order.")
                                ),
                                column(6,

                                       h4("1- Get Sequences") ,
                                       tags$li(" MirCure takes the information about the miRNA sequence from genome (using the provided gff3)."),
                                       tags$li(" Extends the precursor # bases in each extreme in order to allow a correct folding.")

                                )
                              ),
                              fluidRow(
                                column(6,
                                       h4("2- Adjust the annotation"),
                                       tags$li(" MirCure reads the .bam file and adjusts the annotation according to the expression evidence."),
                                       tags$li(" MirCure also saves the previous annotation and compare the final score with the adjust one."),
                                       tags$li("  MirCure would report the structure with a high final score."),
                                       img(src="/appfigs/adjust_annotation.png" ,width="350" , height="300")

                                ),
                                column(6,
                                       h4("3- Fold sequences"),
                                       tags$li(" MirCure calculates the secondary structure of the extended precursor"),
                                       tags$li(" Colors the mature miRNA (red) and star (blue)."),
                                       tags$li(" miRNAs should display a hairpin structure"),
                                       tags$li(" Mature & star secuences should display a 2nts overhang on both extremes."),
                                       h6("~ Score per each extreme: 5 pts= perfect / 3pts = semi-good (3nts overhang)"),
                                       img(src="/appfigs/Bge-Mir-1_pre_fold.jpg" ,width="350" , height="300")
                                ),
                                column(6,
                                       h4("4- Calculate expression"),
                                       tags$li(" After uploading a BAM file with small RNA-seq reads against a genome, caclulate expression of mature/star/prec"),
                                       tags$li(" Creates a bar plot in which for each nt of the precursor shows the number of mapped reads"),
                                       h6("~ Score: > 2 reads in mature = +3 pts // >2 reads in star = +3pts")
                                )
                              ),
                              fluidRow(

                                column(6,
                                       h4("5- Conservation"),
                                       tags$li(" MirCure takes the mature seqeucne and aligns against miRBase"),
                                       tags$li(" Report the alignments (max. 15)"),
                                       h6("~ Score depends on the number of hits. Max=3 "),
                                       img(src="/appfigs/alignments.png" ,width="400" , height="350")

                                ),
                                column(6,
                                       h4("6- Integrate"),
                                       tags$li(" Integrates all the info in a single table"),
                                       tags$li(" Select all candidates you want to accept as real miRNAs (pre-selected > 8pts)"),
                                       tags$li(" Download selected miRNAs")
                                )
                              ),

                              fluidRow(
                                column(6,
                                       h4("7- Give us feedback"),
                                       h5("guillemyllabou-at-gmail-dot-com"),
                                       h5("tianyuan.liu-at-ufl-dot-edu")
                                ),
                                column(6,
                                       h4("Animal / Plant "),
                                       tags$li(" Modifies score system"),
                                       tags$li(" Animal loop length should be shorter than plants"),
                                       tags$li(" In plants it is important to find it conserved, otherwise could be a siRNA")

                                )
                              )
                            )
                          )
                 ),





                 tabPanel(title = 'Sequence Info',
                          # Horizontal line ----
                          tags$hr(),
                          h1("Nucelotide sequences"),

                          h3("miRNA SEQS"),
                          DT::dataTableOutput("mirnaSeqs")
                 ),

                 tabPanel(title = 'RNA Folding',
                          # Horizontal line ----
                          tags$hr(),
                          h1("Secondary structures"),

                          DT::dataTableOutput("mirnaSeqswithplots")

                 ),



                 tabPanel(title = 'Expression Plots',
                          # Horizontal line ----
                          tags$hr(),
                          h1("smallRNA-seq counts"),

                          DT::dataTableOutput("PLOTS"),

                          fluidRow(title = "Expression Plots",
                                   uiOutput("plotouput")
                          )

                 ),


                 tabPanel(title = 'Homology',
                          # Horizontal line ----
                          tags$hr(),
                          h1("Alignments"),

                          DT::dataTableOutput('alignmentsoutput'),
                          hr(),
                          h4("Alignments:"),
                          #tags$style(type='text/css', '#alignments2 {background-color: rgba(255,255,0,0.40); color: green;}'),
                          fluidRow(title = "Aligntab",column(8, align="center",verbatimTextOutput('alignments2')))
                 ),
                 tabPanel(title = 'Integration',
                          # Horizontal line ----
                          tags$hr(),
                          h1("Alignments"),

                          DT::dataTableOutput('Intergartiontable'),
                          hr(),
                          fluidRow(
                            column(width = 12, tableOutput("wsf"),
                                   title = "Aligntab",
                                   column(width =3,
                                          h4("Alignments"),
                                          verbatimTextOutput("Intergartiontable4")
                                   ),
                                   column(width =3,
                                          h4("Precursor Folding"),
                                          htmlOutput("Intergartiontable2")
                                   ),
                                   column(width =6,
                                          h4("Expression Plot"),
                                          uiOutput("Intergartiontable3")
                                   )
                            )
                          ),
                          downloadButton("downloadData", "Download Selected (.csv)"),
                          downloadButton("downloadMatures", "Download Matures (.fa)"),
                          downloadButton("downloadALL", "Download ALL data (csv)"),
                          actionButton("report", "Download the report (pdf)*"),
                          p("*Requires pdflatex")
                 ),
                 tabPanel(title = 'Score stats',
                          # Horizontal line ----
                          tags$hr(),
                          h1("Score Stats"),
                          fluidRow(column(8,align="center",
                                          plotOutput('Scoreshistogram',  width = "80%"))),

                          fluidRow(column(8,align="center",
                                          plotOutput('ScoreBoxplots1',  width = "60%"))),

                          fluidRow(column(8,align="center",
                                          plotOutput('ScoreBoxplots2',  width = "60%")))


                 )
)

