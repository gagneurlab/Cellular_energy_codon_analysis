# This file loads all the csv/txt/xlsx/RDS/etc files needed and displays which
# files have been created/loaded in the workspace.

# Clear enviroment and keep only functions and global/project variables
rm(list = grep(paste(c("^project.variables",
                       lsf.str()),
                     collapse = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Return the memory to the OS
invisible(gc(verbose = FALSE,
             reset = TRUE))

# Create an empty list for all global variables
global.variables <- list()

# Get the command line arguments if they exist
command.line.arguments <- project.variables[["command.line.arguments"]]

# Get a list with all global variables
global.variables <- read.analysis.parameters(command.line.arguments)

# Set the number of threads as a number
global.variables[["threads"]] <- as.numeric(global.variables[["threads"]])

# Set the get.protein.coding amd is.paired.end flags as a logical
global.variables[["get.protein.coding"]] <- as.logical(global.variables[["get.protein.coding"]])

# Set the number of threads as a number
global.variables[["is.paired.end"]] <- as.logical(global.variables[["is.paired.end"]])

# Get the analysis path
analysis.name <- global.variables[["analysis.name"]]

# Get the major condition column
condition.column.major <- global.variables[["condition.column.major"]]

# Get the minor condition column
condition.column.minor <- global.variables[["condition.column.minor"]]

# If the minor condition column is empty, set it to NULL
if(condition.column.minor == "") {
  condition.column.minor <- NULL
  
  global.variables[["condition.column.minor"]] <- NULL
}

# Get the number of threads
number.of.threads <- global.variables[["threads"]]

# Get the samples file name
samples.file.name <- global.variables[["samples.file.name"]]

message("Setting threads to be used to: ", number.of.threads)

# Read the sample annotation file
annotation.table <- fread(global.variables[["annotation.table.path"]])

# Get the number of samples per condition
number.of.samples <- global.variables[["samples.per.condition"]]

# Get the flag for protein coding or not
get.protein.coding <- global.variables[["get.protein.coding"]]

# Get the organism
organism <- global.variables[["organism"]]

# Get the genome version
genome.version <- global.variables[["genome.version"]]

# Get the gencode version
genes.annotation.version <- global.variables[["gencode.version"]]

# Get the conditions of interest
if(global.variables[["conditions.major"]] == "all") {
  
  # If all take all the major conditions
  conditions.major <- annotation.table[, unique(get(global.variables[["condition.column.major"]]))]
  
} else {
  
  # Else split the major conditions between commas
  conditions.major <- unlist(strsplit(global.variables[["conditions.major"]], "[ ]*,[ ]*"))
  
}

# If there is a second condition column as well, get those as well
# Get the conditions of interest
if(is.null(global.variables[["condition.column.minor"]]) == FALSE) {
  if(global.variables[["conditions.minor"]] == "all") {
    
    # If all take all the minor conditions
    conditions.minor <- annotation.table[, unique(get(global.variables[["condition.column.minor"]]))]
    
  } else {
    
    # Else split the minor conditions between commas
    conditions.minor <- unlist(strsplit(global.variables[["condition.column.minor"]], "[ ]*,[ ]*"))
  }
} else {
  conditions.minor <- NULL
}

# Update the major condition  
global.variables[["condition.column.major"]] <- condition.column.major

# Update the minor condition  
global.variables[["condition.column.minor"]] <- condition.column.minor

# Randomly sample N samples from annotation file for each major/minor condition
samples.table <- get.sample.table.from.annotation(annotation.table,
                                                  analysis.name,
                                                  condition.column.major,
                                                  conditions.major,
                                                  condition.column.minor,
                                                  conditions.minor,
                                                  number.of.samples,
                                                  samples.file.name)

# Read the RData file for the Taxonomy object
TxDb <- get.TxDb(organism,
                 genome.version,
                 genes.annotation.version,
                 # overwrite = TRUE,
                 discard.transcripts.types = c("pseudogene", "antisense", "transcribed",
                                               "unprocessed", "retained"),
                 keep.standard = TRUE)

# Get the exons granges
exons.granges <- get.exons.per.gene(TxDb = TxDb,
                                    organism = organism,
                                    genome.version = genome.version,
                                    genes.annotation.version = genes.annotation.version,
                                    # overwrite = TRUE,
                                    protein.coding.only = get.protein.coding)

# Get the introns granges
introns.granges <- get.introns.per.gene(TxDb = TxDb,
                                        exons.by.genes = exons.granges, 
                                        organism = organism,
                                        genome.version = genome.version,
                                        genes.annotation.version = genes.annotation.version, 
                                        protein.coding.only = get.protein.coding,
                                        # overwrite = TRUE,
                                        number.of.threads = number.of.threads)


# Get the exon granges for RNA-seq setup
raw.exons.granges <- get.raw.exons.per.gene(TxDb = TxDb,
                                    organism = organism,
                                    genome.version = genome.version,
                                    genes.annotation.version = genes.annotation.version,
                                    # overwrite = TRUE,
                                    protein.coding.only = get.protein.coding)

# Add the TxDb in the global variables
global.variables[["TxDb"]] <- TxDb

# Add the samples table in the global variables
global.variables[["samples.table"]] <- samples.table

# Add the exon granges in the global variables
global.variables[["exons.granges"]] <- exons.granges

# Add the intron granges in the global variables
global.variables[["introns.granges"]] <- introns.granges

# Add the raw exon granges in the global variables
global.variables[["raw.exons.granges"]] <- raw.exons.granges

message("========== End of load_data.R ==========\n")
