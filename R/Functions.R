runMirCure <- function() {

  ############################
  #### check if Vienna can run
      library("LncFinder")
      if ( class(try(run_RNAfold("AT", RNAfold.path = "RNAfold"),silent = TRUE))== "try-error" ) {
         ViennaDir <- system.file("ViennaRNA",package = "MirCureApp")
        if(.Platform$OS.type == "unix") {
          system(paste0("cd ", ViennaDir," && tar xvzf ViennaRNA-2.4.11.tar.gz" ))
          stop(paste0("For running all the options of MirCure you need to install the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/index.html#download)
To install it UNIX systems, copy and paste the next line in a terminal: \n
cd ", ViennaDir, "/ViennaRNA-2.4.11"," ;  autoconf ; ./configure ; make ; sudo make install"   ))
       } else {# windows
         stop(paste0("For running all the options of MirCure you need to install the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/index.html#download)
For Windows machines, run the .exe file that you will find at: \n", print( ViennaDir) ))
       }
    }else{
      print("RNAfold is present")
      }
  #############################



  appDir <- system.file("shinyAPPs", "MirCure",package = "MirCure")
  if (appDir == "") {
    stop("Could not find shinyAPPs/MirCure directory. Try re-installing `MirCure`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}




installViennaRNA <- function() {
  if(.Platform$OS.type == "unix") {
    #print("unix")

    #ViennaDir <- system.file("ViennaRNA", "ViennaRNA-2.4.11.tar.gz",package = "MirCureApp")
    ViennaDir <- system.file("ViennaRNA",package = "MirCureApp")

  if (ViennaDir == "") {
      stop("Could not find example directory. Try re-installing `MirCureApp`.", call. = FALSE)
    }else{
      #system(paste0(" tar xvzf ", ViennaDir, "/ViennaRNA-2.4.11.tar.gz | ./configure ",ViennaDir," /ViennaRNA-2.4.11 |",  "make ",ViennaDir," /ViennaRNA-2.4.11" ) )
      #system(paste0("cd ", ViennaDir," && tar xvzf ViennaRNA-2.4.11.tar.gz" ))
      #system(paste0( "cd ", ViennaDir,"/ViennaRNA-2.4.11"," && ./configure && make "   ))
      if(file.exists('/usr/local/bin/RNAfold')){"ViennaRNA is present"
        }else{
      system(paste0("cd ", ViennaDir," && tar xvzf ViennaRNA-2.4.11.tar.gz" ))
      print("For running all the options of MirCure you need to install the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/index.html#download)",  quote = FALSE)
      print("For UNIX systems, copy and paste the next line in a terminal: ",  quote = FALSE)

      print(paste0( "cd ", ViennaDir, "/ViennaRNA-2.4.11"," && ./configure && make && sudo make install"   ))
      #print(paste0( "cd ", ViennaDir, " && tar xvzf ViennaRNA-2.4.11.tar.gz && cd ViennaRNA-2.4.11"," && ./configure && make && sudo make install"   ))
        }
      }
    } else {# windows
      print("For running all the options of MirCure you need to install the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/index.html#download)",  quote = FALSE)
      print("For Windows machines, run the .exe file that you willl find at: ",  quote = FALSE)
      print( ViennaDir)
  }
}
