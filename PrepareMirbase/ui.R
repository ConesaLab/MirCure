# Define UI for dataset viewer app ----
library("shinythemes")

ui <- navbarPage("Prepare miRBase for MirCure",theme = shinytheme("sandstone"),

tabPanel(title='Load Data',
      # Sidebar layout with a input and output definitions ----
      sidebarLayout(
        #position = "right",
        # Sidebar panel for inputs ----
        sidebarPanel(

          fileInput("inputfile1", "Mirbase (.gff3)",
                    multiple = FALSE,
                    accept = c(".gff3")),

        actionButton("run", "RUN!"),
        hr(),
        downloadButton("downloadprecs", "Download Precursors Gff3"),
        hr(),
        downloadButton("downloadp5p", "Download 5' arms Gff3"),
        hr(),
        downloadButton("download3p", "Download 3' arms Gff3")


           ),

# Main panel for displaying outputs ----
      mainPanel(
        h1("Hello world!"),
        h3("Input data"),
        tableOutput("visualize"),
        h3("Precursors"),
        tableOutput("precs"),
        h3("5' arms miRNAs"),
        tableOutput("fiveP"),
        h3("3'arms miRNAs"),
        tableOutput("threeP")

)
)
)
)
