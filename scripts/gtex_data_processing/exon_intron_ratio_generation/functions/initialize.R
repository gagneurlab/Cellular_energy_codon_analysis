# This file loads all the packages,libraries and data needed regarding the
# workspace, loads the functions.R script and sets the global variables.

# Jump-start analysis for gtex
# project.variables[["command.line.arguments"]] <- c("--analysis.name", "all-samples-all-tissues")

# Check that the command line arguments are valid
valid.arguments <- check.command.line.arguments(project.variables[["command.line.arguments"]])

# Print a message for analysis begin
message("Start mrna half-life generation pipeline at: ", date())

# Set seed for this session
set.seed(24101992) 

# Clear enviroment
rm(list = grep(paste("^project.variables",
                     "^check.and.install.packages",
                     "^check.command.line.arguments",
                     "^check.parameters.validity",
                     sep = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Return the memory to the OS
invisible(gc(verbose = FALSE,
             reset = TRUE)) 

# Project Packages to be installed
cran.packages <- c( "data.table",
                    "here",
                    "R.utils",
                    "magrittr",
                    "stringr",
                    "pheatmap",
                    "RColorBrewer",
                    "parallel",
                    "ggplot2",
                    "ggrepel",
                    "ggpp",
                    "ggpubr",
                    "gtable",
                    "grid",
                    "ggplotify",
                    "gtools",
                    "ggsignif",
                    "glmnet",
                    "splitstackshape",
                    "plotly",
                    "parallelDist",
                    "yaml",
                    "retry",
                    "pbmcapply",
                    # "ProliferativeIndex",
                    "tibble")

bioconductor.packages <- c( "rtracklayer",
                            "GenomicAlignments",
                            "GenomicRanges",
                            "GenomicFeatures",
                            "Rsamtools",
                            "biomaRt",
                            "DESeq2",
                            "IRanges",
                            "BiocParallel",
                            "BiocIO",
                            "fgsea",
                            "EnrichmentBrowser",
                            "org.Hs.eg.db",
                            "org.Mm.eg.db",
                            "org.At.tair.db",
                            "GO.db",
                            "coRdon",
                            "ViSEAGO",
                            "limma",
                            "BSgenome.Athaliana.TAIR.TAIR9",
                            "BSgenome.Hsapiens.UCSC.hg19",
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "BSgenome.Mmusculus.UCSC.mm9",
                            "BSgenome.Mmusculus.UCSC.mm10",
                            "motifStack")

github.packages <- c("flashpcaR", "fastglm")

# Add packages used during development only
if(project.variables[["development.stage"]] == TRUE & project.variables[["use.packrat"]] == TRUE) {
  cran.packages <- c(cran.packages, "rbenchmark")
  
  # Temporary reset current working directory
  # in order to work packrat package installation
  setwd(here())
  
  # Install missing packages
  check.and.install.packages(cran.packages, origin = "CRAN")
  
  # Install missing packages
  check.and.install.packages(bioconductor.packages, origin = "Bioconductor")
  
  # Install missing packages
  check.and.install.packages(github.packages, origin = "GitHub")
  
  # Reset current working
  setwd(here("src"))
}

# Combine all packages regardless of origin in a vector
packages <- c(cran.packages, bioconductor.packages, github.packages)

message("Loading packages...")

# Load all packages
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

# Save loaded packages in packrat
if (project.variables[["development.stage"]] == TRUE & project.variables[["use.packrat"]] == TRUE) {
  
  cat("Check Packrat packages status and snapshot if needed...\n")
  
  # Temporary reset current working directory
  # in order to work packrat package installation
  setwd(here())
  
  # Check the packrat packages status
  packages.status <- status()
  
  # Are the packages in the last snapshot, the same as in the local downloaded packages?
  packages.are.up.to.date <- all(packages.status$packrat.version == packages.status$library.version)
  
  # If not, snapshot the loaded packages
  if (packages.are.up.to.date == FALSE | is.na(packages.are.up.to.date) == TRUE) {
    snapshot()  
  }
  
  # Reset current working
  setwd(here("src"))
}

source("functions.R")

message("========== End of initialize.R ==========\n")
