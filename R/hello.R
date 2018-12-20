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
