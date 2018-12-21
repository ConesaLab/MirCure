# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



##################https://deanattali.com/2015/04/21/r-package-shiny-app/
hello <- function() {
  print("Hello, world!")
}



runExample <- function() {
  appDir <- system.file("shinyAPPs", "basicshinnyapp",package = "mirQCApp")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mirQCApp`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


runMirPlot <- function() {
  appDir <- system.file("shinyAPPs", "mirPlot",package = "mirQCApp")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mirQCApp`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}




installViennaRNA <- function() {
  if(.Platform$OS.type == "unix") {
    #print("unix")

    #ViennaDir <- system.file("ViennaRNA", "ViennaRNA-2.4.11.tar.gz",package = "mirQCApp")
    ViennaDir <- system.file("ViennaRNA",package = "mirQCApp")

  if (appDir == "") {
      stop("Could not find example directory. Try re-installing `mirQCApp`.", call. = FALSE)
    }else{
      #system(paste0(" tar xvzf ", ViennaDir, "/ViennaRNA-2.4.11.tar.gz | ./configure ",ViennaDir," /ViennaRNA-2.4.11 |",  "make ",ViennaDir," /ViennaRNA-2.4.11" ) )
      #system(paste0("cd ", ViennaDir," && tar xvzf ViennaRNA-2.4.11.tar.gz" ))
      #system(paste0( "cd ", ViennaDir,"/ViennaRNA-2.4.11"," && ./configure && make "   ))
      print("For running all the options of MirQC you need to isntall the ViennaRNA",  quote = FALSE)
      print("For UNIX systems, copy and paste the next line in a terminal: ",  quote = FALSE)
      print(paste0( "cd ", ViennaDir, " && tar xvzf ViennaRNA-2.4.11.tar.gz && cd ViennaRNA-2.4.11"," && ./configure && make && sudo make install"   ))
    }
    } else {# windows
      print("For running all the options of MirQC you need to isntall the ViennaRNA",  quote = FALSE)
      print("For Windows machines, run the .exe file that you willl find at: ",  quote = FALSE)
      print( ViennaDir)
  }
}
