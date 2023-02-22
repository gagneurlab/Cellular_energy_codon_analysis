# Sets the current working directory and calls the following scripts:
# initialize.R, load_data.R, pull_data_from_DB.R, build.R, analyze.R

# Clear enviroment
rm(list = ls())

# Return the memory to the OS
invisible(gc(verbose = FALSE,
             reset = TRUE))

# Am I on a development stage?
development.stage <- TRUE

# Should I use packrat?
use.packrat <- FALSE

# If yes, should i reset it?
reset.packrat <-  TRUE

# Try to get the analysis parameters from the command line
arguments <- commandArgs(trailingOnly = TRUE)

if(is.null(arguments) == FALSE) {
  
  # Set the flag to use arguments from the command line
  use.command.line.parameters <- TRUE
  
} else {
  # If no command line arguments were provided, try to use a file later 
  use.command.line.parameters <- FALSE
}

# Set project variables
project.variables <- list("development.stage" = development.stage,
                          "use.packrat" = use.packrat,
                          "use.command.line.parameters" = use.command.line.parameters,
                          "command.line.arguments" = arguments)

# Project Packages to be installed
cran.packages <- c("here", "devtools")
github.packages <- c("dokato/todor")

# Source the utilities file
source("src/functions_utilities.R")

# Should I use the packrat enviroment or my personal/lab package environment 
if(use.packrat == TRUE) {
  
  # If yes, does the packrat folder exist or we should create it?
  if ((dir.exists("packrat") == FALSE | reset.packrat == TRUE) & 
      project.variables[["development.stage"]] == TRUE) {
    
    # Set mirrors and repositories
    chooseCRANmirror(graphics = FALSE, ind = 1)

    chooseBioCmirror(graphics = FALSE, ind = 1)

    setRepositories(graphics = FALSE, ind = 1:4)
    # Remove old packrat
    system("rm -rf packrat/ .Rprofile")
    
    # Install packrat as the base package manager
    check.and.install.packages("packrat")
    
    # Load Packrat
    library(packrat)
    
    # Initialize packarat
    init(getwd(), restart = FALSE)
    
    ######################################
    # After init it restarts rsession
    # and nothing after here is executed
    ######################################

    # Set packrat mode ON
    packrat_mode(on = TRUE)
    
    # Install the Project CRAN Packages needed for development
    check.and.install.packages(cran.packages)
    
    # Load the project CRAN packages
    invisible(lapply(cran.packages, require, character.only = TRUE))
    
    # Load github packages
    check.and.install.packages(github.packages, origin = "GitHub")
    
    # Take a packrat snapshot
    snapshot()
    
  } else {
    #if(use.packrat == TRUE) {
      library("packrat")

      # # Set packrat mode ON
      packrat_mode(on = TRUE)
      
      if (all(cran.packages %in% .packages()) == FALSE) {
        # Install the Project CRAN Packages needed for development
        check.and.install.packages(cran.packages)
        
        # Load all the packages (packrat, here, devtools)
        suppressPackageStartupMessages(invisible(lapply(cran.packages, library, character.only = TRUE)))
        
        # Load github packages
        check.and.install.packages(github.packages, origin = "GitHub")
      }
      
      # Load all the packages (packrat, here, devtools)
      suppressPackageStartupMessages(invisible(lapply(cran.packages,
                                                      library,
                                                      character.only = TRUE)))
  }
} else {

  system("rm -rf .Rprofile .Rhistory")
  
    # If we use the personal/lab package environment, just load the packages
  suppressPackageStartupMessages(invisible(lapply(cran.packages,
                                                  library,
                                                  character.only = TRUE)))
}

# Scripts to call
files.to.load <- c( "initialize.R",
                    "load_data.R",
                    "build.R")#,
                    # "analyze.R")

# Set the currenct working directory
setwd(here("src"))

message("========== End of main.R ==========\n")

# Run the whole analysis
analysis <- lapply(files.to.load, source)

