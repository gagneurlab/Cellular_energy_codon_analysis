# ==== Load functions ====

read.analysis.parameters <- function( command.line.arguments = NULL,
                                      analysis.parameters.file = file.path(here(),
                                                                           "data-input",
                                                                            "analysis-parameters.csv")) {
  #
  # Makes the list of global variables with the analysis parameters
  #
  # Requires:
  #   data.table
  #   here
  #
  # Load command:
  # sapply(c("data.table", "here"), library, character.only = T)
  #
  # Args:
  #   command.line.arguments:   Default is NULL. A vector with the command line parameters
  #   analysis.parameters.file: Default is the "data-input/analysis-parameters.csv".
  #                             The file path to the analysis parameters' csv
  #
  # Returns:
  #   A list with the analysis parameters
  #
  
  # Prepare the analysis parameters list
  analysis.parameters <- list()
  
  # We can either read the arguments for the analysis from the command line
  # or by default from the data-input/analysis-parameters.csv file
  # If nothing is provided, stop the analysis
  
  if((is.null(command.line.arguments) == TRUE | length(command.line.arguments) == 0) &
     file.exists(analysis.parameters.file) == TRUE) {
    
    message("Reading analysis-parameters.csv from data-input/.")
    
    # Read the analysis parameters table
    analysis.parameters.table <- fread(analysis.parameters.file)
    
    # Put the argument values on a list
    analysis.parameters <- as.list(analysis.parameters.table[, value])
    
    # Put the parameter names as names on the analysis parameters list
    names(analysis.parameters) <- analysis.parameters.table[, parameter]
    
    # Check the parameter validity
    check.parameters.validity(analysis.parameters)
    
  } else if(is.null(command.line.arguments) == FALSE  & 
            length(command.line.arguments) > 0) {
    
    # Check that i have valid arguments
    analysis.parameters <- check.command.line.arguments(arguments = command.line.arguments,
                                                         tolerate.no.arguments = FALSE)
    
  }
  
  # Check if the directory already exists or if it is a novel run
  analysis.directory <- file.path(here(),
                                  "data-input",
                                   analysis.parameters[["analysis.name"]])
  
  # Depending on if one or more arguments are provided,
  # Check that the analysis folder exists
  if(length(analysis.parameters) != 1 & 
     dir.exists(analysis.directory) == FALSE) {
    
    # Check the validity of the parameters
    check.parameters.validity(analysis.parameters)
    
    message("mRNA half-life analysis directory ",
            analysis.parameters[["analysis.name"]],
            " does not exist... Lets create it!")
    
    # Create the analysis folder
    dir.created <- dir.create(analysis.directory)
    
    # Store the analysis directory
    analysis.parameters[["analysis.input.directory"]] <- analysis.directory
    
    # Create the new parameter table
    new.analysis.parameters.table <- data.table("parameter" = names(analysis.parameters),
                                                "value" = analysis.parameters)
    
    # Prepare the new path
    new.analysis.parameters.table.path <- file.path(analysis.directory,
                                                    "analysis-parameters.csv")
    
    # And now write the new table
    fwrite(new.analysis.parameters.table,
           new.analysis.parameters.table.path,
           sep = ",")
    
    return(analysis.parameters)
    
  } else if(length(analysis.parameters) == 1 & 
            dir.exists(analysis.directory) == FALSE) {
    
    stop("--analysis.name folder does not exist.\n",
         "Please specify all the arguments through the command line or through a file.",
         "Available parameters are the following:\n",
         paste0("\t", "--", check.command.line.arguments(report.arguments.only = TRUE),"\n"),
         call. = FALSE)
    
  } else {
    
    # Prepare the path
    analysis.parameters.file <- file.path(analysis.directory,
                                          "analysis-parameters.csv")
    
    message("Scanning for ", analysis.directory,"...\n",
            "File already exists. Reading from analysis directory...")
    
    # Read the analysis parameters table
    analysis.parameters.table <- fread(analysis.parameters.file)
    
    # Put the argument values on a list
    analysis.parameters <- as.list(analysis.parameters.table[, value])
    
    # Put the parameter names as names on the analysis parameters list
    names(analysis.parameters) <- analysis.parameters.table[, parameter]
    
    # Check the validity of the parameters
    check.parameters.validity(analysis.parameters)
    
    return(analysis.parameters)
    
  }
}

get.sample.table.from.annotation <- function(annotation.table,
                                            analysis.name,
                                            condition.column.major, 
                                            conditions.major,
                                            condition.column.minor = NULL, 
                                            conditions.minor = NULL,
                                            number.of.samples = 5,
                                            samples.file.name = "samples-to-run.csv",
                                            overwrite = FALSE) {
  #
  # Reads the annotation table of GTEX containing information for all the experiments,
  # and draws N random samples for the tissues/subtissues of interest
  #
  # Args:
  #   annotation.table:       A data.table containing the metadata of GTEX samples
  #   analysis.path           The analysis path string
  #   condition.column.major: The column name of the major condition
  #   conditions.major:       The major conditions of interest
  #   condition.column.minor: Default is NULL. The column name of the minor condition
  #   conditions.minor:       Default is NULL. The minor conditions of interest
  #   number.of.samples:      Default is 5. The number of random samples to draw or "all" to get 
  #                           all samples from a tissue
  #   samples.file.name:      Default is "gtex-samples.csv". The file name with the
  #                           gtex random samples
  #   overwrite:              Default is FALSE. Should I overwrite the sample annotation table?
  #
  # Returns:
  #   A data.table containing N random samples of each tissue/subtissue
  #
  
  # Prepare the path of the output file
  samples.file.path <- file.path(here(),
                             "data-input",
                             analysis.name,
                                  samples.file.name)
  
  if(file.exists(samples.file.path) == FALSE | overwrite == TRUE) {
    
    message(samples.file.path," not found, it will be generated...\n")
    
    # Now from all the rows of desired conditions, sample N for each combination
    if(number.of.samples == "all"){
      
      # Do I have only a major condition or a minor as well.
      # If I have 2, get samples for both
      if(is.null(condition.column.minor) == TRUE) {
        samples.ids <- annotation.table[  get(condition.column.major) %in% conditions.major,
                                          RNA_ID,
                                          by = condition.column.major][, RNA_ID]
        
      } else {
        samples.ids <- annotation.table[  get(condition.column.major) %in% conditions.major &
                                            get(condition.column.minor) %in% conditions.minor,
                                          RNA_ID,
                                          by = c(condition.column.major, condition.column.minor)][, RNA_ID]
      }
      
    } else {
      
      # Do I have only one condition of interest column, or 2.
      # If I have 2, get samples for both
      if(is.null(condition.column.minor) == TRUE) {
        samples.ids <- annotation.table[  get(condition.column.major) %in% conditions.major,
                                          .(RNA_ID = sample(RNA_ID, replace = FALSE, size = number.of.samples)),
                                          by = condition.column.major][, RNA_ID]
        
      } else {
        samples.ids <- annotation.table[  get(condition.column.major) %in% conditions.major &
                                            get(condition.column.minor) %in% conditions.minor,
                                          .(RNA_ID = sample(RNA_ID, replace = FALSE, size = number.of.samples)),
                                          by = c(condition.column.major, condition.column.minor)][, RNA_ID]
      }
    
    }
    
    # Add a new column to the annotation table with the bam files to be sampled
    annotation.table[, is.sampled :=  RNA_ID %in% samples.ids]
    
    # Does SRA_Sample column exist? If not, create a dummy one
    if(!"SRA_Sample" %in% colnames(annotation.table)) {
      annotation.table[, SRA_Sample := ""]
    }
    
    # Does sex column exist? If not, create a dummy one
    if(!"sex" %in% colnames(annotation.table)) {
      annotation.table[, sex := ""]
    }
    
    # Prepare the sample table for export
    samples.table <- annotation.table[is.sampled == TRUE, .SD, 
                                                          .SDcols = c("RNA_ID",
                                                            "SRA_Sample",
                                                            condition.column.major,
                                                            condition.column.minor,
                                                            "RNA_BAM_FILE",
                                                            "sex")]
    
    message("Write the ", samples.file.path," file...\n")
    
    # Write the table with the SRA ids of interest
    fwrite(samples.table,
           samples.file.path,
           quote = FALSE,
           sep = ",")
    
  } else {
    
    message("Read the samples csv file.\n")
    
    # Read the gtex samples file
    samples.table <- fread(samples.file.path)
  }
  
  return (samples.table)
}

get.gencode.annotation <- function(gencode.organism = "human",
                                   gencode.version = "gencode.v19") {
  #
  # Downloads the gene annotation for a specific version
  #
  # Requires:
  #
  # Load command:
  # 
  #
  # Args:
  #   organism:         Default is "Homo Sapiens". The gencode organism of interest.
  #                     Can be Homo Sapiens or Mus Musculus    
  #   gencode.version:  Default is "gencode.v19". The gencode version to download e.g.
#                       "gencode.v38"
  #
  # Returns:
  #   The file path of the gene annotation
  #
  
  # Set the file name of the gene annotation as found on the site of genecode
  gencode.annotation.file <- paste0(gencode.version,
                                    ".annotation.gtf.gz")
  
  # Make the full path where it will be stored
  gencode.annotation.file.path <- file.path(here(),
                                        "data-input",
                                        "gencode-annotations",
                                        paste0(gencode.organism, ".",
                                               gencode.annotation.file))
  
  # Get the gencode version number
  gencode.version.number <- unlist(strsplit(gencode.version, ".v"))[2]
    
  # Prepare the ftp path
  ftp.path <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_", gencode.organism,"/release_",
                     gencode.version.number,
                     "/",
                     gencode.annotation.file)
  
  # If the file does not exist on the data-output folder, download it there
  if(file.exists(gencode.annotation.file.path) == FALSE) {
    download.file(ftp.path,
                  destfile = gencode.annotation.file.path)
    
  }
  
  # Return the downloaded file
  return (gencode.annotation.file.path)
  
}

get.TxDb <- function( organism,
                      genome.version,
                      genes.annotation.version,
                      discard.transcripts.types = c("pseudogene", "antisense", "transcribed",
                                                      "unprocessed", "retained"),
                      keep.standard = TRUE,
                      overwrite = FALSE) {
  #
  # Reads and returns the sqlite file for the TxDb object for genome version X gencode version Y.
  # If this file does not exists, it creates it and then returns it.
  #
  # Requires:
  #   data.table
  #   here
  #   GenomicFeatures
  #   GenomicRanges
  #   BiocIO
  #
  # Load command:
  #   sapply(c("data.table", "here", "GenomicFeatures", "GenomicRanges", "BiocIO"), 
  #           library, character.only = T)
  #
  # Args:
  #   organism:                 The name of the organism of interest e.g. "Homo sapies"
  #                             or "Arabidopsis thaliana"
  #   genome.version:           The genome version to use e.g. gencode.v38
  #   genes.annotation.version: The gene annotation version to use
  #   discard.transcripts.types:Default is c("pseudogene", "antisense", "transcribed",
  #                             "unprocessed", "retained"). The transcript types to discard
  #   keep.standard:            Default TRUE. Keep only standard chromosomes or not
  #   overwrite:                Default is FALSE. Should I overwrite the TxDb object?
  #
  # Returns:
  #   A TxDb object for genome of interest and the gencode version of interest
  #
  
  if(organism == "Homo sapiens" |
     organism == "Mus musculus") {
    
    if(organism == "Homo sapiens") {
      
      gencode.organism <- "human"
    } else if(organism == "Mus musculus") {
      gencode.organism <- "mouse"
    }
    
    # Turn the organism to lower case
    organism.lower <- gsub(" ", "\\.", tolower(organism))
    
    # Get the gencode annotation path for the version of interest
    genes.annotation.path <- get.gencode.annotation(gencode.organism,
                                                    genes.annotation.version)
    
    # Prepare the TxDb database name
    TxDb.sqlite.name <- paste("TxDb",
                              organism.lower,
                              genome.version,
                              genes.annotation.version,
                              "sqlite",
                              sep = ".")
    
  } else {
    
    # Turn the organism to lower case
    organism.lower <- gsub(" ", "\\.", tolower(organism))
    
    # Get the gencode annotation path for the version of interest
    genes.annotation.path <- file.path(here(),
                                   "data-input",
                                   "arabidopsis-gene-annotations",
                                   "Araport11_GFF3_genes_transposons.Mar92021.cleanedFor.makeTxDbFromGFF.gff")
    
    # Prepare the TxDb database name
    TxDb.sqlite.name <- paste("TxDb",
                              organism.lower,
                              genome.version,
                              genes.annotation.version,
                              "sqlite",
                              sep = ".")
  }
  
  # Create the output path for the TxDb
  TxDb.sqlite.path <- file.path(here(),
                            "data-input",
                            "TxDb-objects",
                            TxDb.sqlite.name)
  
  # If the TxDb object does not exists make it 
  if(file.exists(TxDb.sqlite.path) == FALSE | overwrite == TRUE) {
    
    message("No TxDb object found...")
    message("Create the TxDb object for ", organism,
            ", genome version ",
            genome.version,
            ", gene annotation version ",
            genes.annotation.version,".\n") 
    
    # Make a TxDb object from the gtf file
    # TxDb <- makeTxDbFromGFF(genes.annotation.path,
    #                         organism = organism)
    
    message("Read the gtf file as a GRanges object...")
    
    # Read the gtf as Granges
    gtf.granges <- import(genes.annotation.path)
    
    if(is.null(discard.transcripts.types) == FALSE) {
      
      message("Discard transcripts of the following types: ",
              paste(discard.transcripts.types, collapse = ", "), "...")
      
      # Prepare the pattern
      transcript.pattern <- paste("(", discard.transcripts.types,")",
                                  sep = "", collapse = "|")
      
      # Discard transcripts belonging to one bad category type
      transcripts.to.keep <- grep(transcript.pattern,
                                       gtf.granges$transcript_type,
                                       invert = TRUE)
      
      # Now get the transcripts GRanges for all the relevant transcripts
      gtf.granges.cleaned <- gtf.granges[transcripts.to.keep]
      
      message("Create the TxDb object...")
      
      # And prepare the TxDb object
      TxDb <- makeTxDbFromGRanges(gr = gtf.granges.cleaned)                        
    } else {
      message("Create the TxDb object...")
      # Prepare the TxDb object
      TxDb <- makeTxDbFromGRanges(gr = gtf.granges) 
    }
    
    if(keep.standard == TRUE) {
      
      message("Keep only standard chromosomes...\n")
      
      # Keep only the meaningfull chromosomes
      TxDb.standard.chromosomes <- keepStandardChromosomes(TxDb,
                                                           species = gsub(" ",
                                                                          "_",
                                                                          organism))
      
      # Update the TxDb object
      TxDb <- TxDb.standard.chromosomes
      
    } else {
      
      message("Keep all chromosomes...\n")
    } 
    
    # And save them in a sqlite file for faster analysis
    saveDb(TxDb, 
           file =  TxDb.sqlite.path)
  } else {
    message("TxDb object found!\n")
    message("Loading TxDb object for ", organism,
            ", genome version ",
            genome.version,
            ", gene annotation version ",
            genes.annotation.version,".\n")
    
    # If the file already exists read it directly from the RDdata file
    TxDb <- loadDb(TxDb.sqlite.path)
  }
  
  # Return the TxDb object
  return (TxDb)
}

get.raw.exons.per.gene <- function(TxDb,
                                   organism,
                                   genome.version,
                                   genes.annotation.version,
                                   protein.coding.only = FALSE, 
                                   overwrite = FALSE) {
  #
  # Loads a GRangesList object with the raw exons grouped by gene and the
  #   names of the genes if it exists, otherwise it creates it
  #
  # Requires:
  #   data.table
  #   here
  #   GenomicFeatures
  #   biomaRt
  #   magrittr
  #   GenomicRanges
  #
  # Load command:
  # sapply(c("data.table", "here", "GenomicFeatures", "biomaRt", "magrittr",
  #       "GenomicRanges"), library, character.only = T)
  #
  # Args:
  #   TxDb:                     The database TxDb object
  #   organism:                 The name of the organism of interest e.g. "Homo sapies" or "Arabidopsis thaliana"
  #   genome.version:           The genome version to use e.g. gencode.v38
  #   genes.annotation.version: The gene annotation version to use
  #   protein.coding.only:      Default is FALSE. Should i get only the coding genes
  #   overwrite:                Default is FALSE. Should I overwrite the file?
  #
  # Returns:
  #   A GRangesList object with the exons grouped by gene and the
  #   names of the genes
  #
  
  # Read the TxDb object for e.g. hg19 v19
  # Beware that chomosomes are stored as 
  # chr1, chr2, ..., chrX, chrY, chrMT
  
  # Prepare an empty total granges object
  raw.exons <- NULL
  
  # Add the gencode prefix for human and mouse
  if(organism == "Homo sapiens" |
     organism == "Mus musculus") {
    
  }
  
  # Prepare the folder path for the total granges file
  raw.exons.folder <- file.path(here(),
                                  "data-input",
                                  "exon-coordinates",
                                  "raw-exons")
  
  # Make the folder if it does not exist
  if(dir.exists(raw.exons.folder) == FALSE) {
    dir.create(raw.exons.folder)
  }
  
  # Prepare the file path for the total granges file
  raw.exons.file <- file.path(raw.exons.folder,
                                paste0(genome.version, "-", genes.annotation.version, "-",
                                       ifelse(protein.coding.only == TRUE,
                                              "raw-exons-protein-coding-only.rds",
                                              "raw-exons.rds")))
  
  # If the exons file does not exist, create it, otherwise read it
  if(file.exists(raw.exons.file) == FALSE | overwrite == TRUE) {
    
    message("Raw exons file not found. Let's create it...\n")
    
    # Group the exons by gene
    exons.by.genes <- exonsBy(TxDb, by = "gene")
    
    # Get the exons by transcript
    exons.by.transcript <- exonsBy(TxDb, by = "tx", use.names= TRUE)
    
    # Read the table with the hosts
    organisms.to.hosts <- fread(here("data-input/organisms-to-hosts.csv"))
    
    # Save again to avoid mixing
    genome.version.input <- genome.version
    
    # Get the hostname
    hostname <- organisms.to.hosts[genome.assembly == genome.version.input, host]
    
    # Get the hostname
    mart.to.access <- organisms.to.hosts[genome.assembly == genome.version.input, mart]
    
    # Prepare to connect to Ensembl
    ensembl <- useMart(mart.to.access, host = hostname)
    
    # Get the dataset of interest
    dataset.of.interest <- listDatasets(ensembl) %>%
      as.data.table %>%
      .[version == organisms.to.hosts[genome.assembly == genome.version.input,
                                      genome.version],
        dataset]
    
    # Use the dataset
    ensembl.mart <-  useDataset(dataset.of.interest, mart = ensembl)
    
    # Get the ensembl gene column
    ensembl.gene.column <- ifelse(genome.version.input == "TAIR10",
                                        "ensembl_gene_id",
                                        "ensembl_gene_id_version")
    
    # Get the gene name according to the organism
    if(organism == "Homo sapiens") {
      gene.name.column <- "hgnc_symbol"
    } else if(organism == "Mus musculus") {
      gene.name.column <- "mgi_symbol"
    } else {
      stop("Peos in raw exons")
    }
    
    # Get the the information for all genes, all chromosomes as well the
    # type of transcript (protein coding or whatever)
    # and the description of the gene
    ensembl.genes <- data.table(getBM(attributes = c(ensembl.gene.column,
                                                     "ensembl_transcript_id_version",
                                                     "chromosome_name",
                                                     "description",
                                                     "transcript_biotype"),
                                      mart = ensembl.mart))
    
    # Make a vector with all the standard chromosomes
    # In Ensembl chromosomes are stored as 1, 2, ..., X, Y, MT
    standard.chromosomes <- standardChromosomes(TxDb, 
                                                species = gsub(" ",
                                                               "_",
                                                               organism))
    
    # Make the standard chromosomes ensemble style
    if(all(grepl("(chr[0-9MCYX]+)", standard.chromosomes, ignore.case = T)) == TRUE) {
      
      # Get the chromosomes Ensembl style
      standard.chromosomes <- gsub("chr", "", standard.chromosomes, ignore.case = T)
      
      # Update the Chloroplast chromosome names
      ensembl.genes[chromosome_name == "Pt", chromosome_name := "C"]
      
      # Update the Mitochondria chromosome names
      ensembl.genes[chromosome_name == "Mt" |
                      chromosome_name == "MT", chromosome_name := "M"]
    }
    
    # Now subset the Ensembl data table by keeping only the protein coding
    # genes, and only those which belong to a standard chromosome
    if(protein.coding.only == TRUE) {
      ensembl.genes <- ensembl.genes[ transcript_biotype == "protein_coding" &
                                        chromosome_name %in% standard.chromosomes, ]
      
      
    } else {
      ensembl.genes <- ensembl.genes[chromosome_name %in% standard.chromosomes, ]
    }
    
    # Get the names of the genes for ensemble
    ensembl.genes[, ensembl_gene_id := gsub("\\.[0-9]+", "",
                                                          get(ensembl.gene.column))]
    
    # Get the gene names
    exons.by.gene.ids <- gsub("\\.[0-9]+", "", names(exons.by.genes))
    
    message("Get the total counts coordinates...")
    
    # Now subset the genes to get only the common between genome annotation/ensembl
    exons.by.genes <- exons.by.genes[exons.by.gene.ids %in%
                                       ensembl.genes[, ensembl_gene_id]]
    
    message("Save the raw exons coordinates as RDS...")
    
    # Save the total granges in a file
    saveRDS(exons.by.genes,
            raw.exons.file)
    
    # And prepare them for export
    raw.exons <- exons.by.genes
    
    message("Raw exons file created!")
  } else {
    message("Raw exons file found! Loading...\n")
    
    # Read the exons file
    raw.exons <- readRDS(raw.exons.file)
  }
  
  return (raw.exons)
}

discard.big.exons.from.granges <- function(flanked.exons.by.genes,
                                           reduced.exons.by.genes,
                                           max.exon.size.percentage = 80,
                                           number.of.threads = 8) {
  
  #
  # Loads a GRangesList object with the exons grouped by gene and the
  #   names of the genes if it exists, otherwise it creates it
  #
  # Requires:
  #   data.table
  #   here
  #   GenomicFeatures
  #   biomaRt
  #   magrittr
  #   GenomicRanges
  #
  # Load command:
  # sapply(c("data.table", "here", "GenomicFeatures", "biomaRt", "magrittr",
  #       "GenomicRanges"), library, character.only = T)
  #
  # Args:
  #   TxDb:                     The database TxDb object
  #   organism:                 The name of the organism of interest e.g. "Homo sapies" or "Arabidopsis thaliana"
  #   genome.version:           The genome version to use
  #   genes.annotation.version: The gene annotation version to use
  #   protein.coding.only:      Default is FALSE Get the exons for the protein coding genes only
  #   overwrite:                Default is FALSE. Should I overwrite the file?
  #
  # Returns:
  #   A GRangesList object with the exons grouped by gene and the
  #   names of the genes
  #
  
  # Read the TxDb object for e.g. hg19 v19
  # Beware that chomosomes are stored as 
  # chr1, chr2, ..., chrX, chrY, chrMT
  
  # Prepare the filtered exons
  exons.by.genes.filtered <- mcmapply(function(exon.grages,
                                               max.exon.size.percentage) {
                                    
                                    # Get the starting point of the granges
                                    grange.start <- head(start(exon.grages), n = 1)
                                    
                                    # And the end point of the granges
                                    grange.end <- tail(end(exon.grages), n = 1)
                                    
                                    # Prepare the granges total width
                                    grange.length <- grange.end - grange.start + 1
                                    
                                    # Add the total length of the exonic ranges 
                                    exon.grages$total.reduced.width = rep(grange.length,
                                                                          length(exon.grages))
                                    
                                    # Add the length of each exon
                                    exon.grages$exon.width <-  unlist(width(exon.grages))
                                    
                                    # Calculate the percentage of each exon in the total exonic length
                                    exon.grages$exon.percentage <-  exon.grages$exon.width * 100 /
                                                                    exon.grages$total.reduced.width
                                    
                                    # Filter out those big exons
                                    exon.grages.cleaned <- exon.grages[exon.grages$exon.percentage <= max.exon.size.percentage]
                                    
                                    exon.grages.big.exons <- exon.grages[exon.grages$exon.percentage > max.exon.size.percentage]
                                    
                                    return (list(list("cleaned" = exon.grages.cleaned,
                                                      "with.big.exons" = exon.grages.big.exons)))
                                    },
                                    exons.by.genes[1:1000],
                                    # sum(width(reduced.exons.by.genes))[28577:28578],
                                    MoreArgs = list(max.exon.size.percentage = max.exon.size.percentage),
                                    mc.cores = number.of.threads)

  # Get the cleaned ranges
  exons.by.genes.cleaned <- mclapply(exons.by.genes.filtered,
                                             function(exon.grages) exon.grages$cleaned)
  
  # Get the genes with discarded exons
  exons.by.genes.with.big.exons <- mclapply(exons.by.genes.filtered,
                                             function(exon.grages) exon.grages$with.big.exons)
  
  # Set the results as a granges list
  exons.by.genes.cleaned <- GRangesList(exons.by.genes.cleaned)
  
  # Set the results as a granges list with big exons
  exons.by.genes.with.big.exons <- GRangesList(exons.by.genes.with.big.exons)
  
  # Get the genes with big exons
  genes.with.big.exons <- names(unlist(width(exons.by.genes.with.big.exons)) > 0)
  
  return (exons.by.genes.cleaned)
}

get.exons.per.gene <- function(TxDb,
                               organism,
                               genome.version,
                               genes.annotation.version,
                               protein.coding.only = FALSE, 
                               overwrite = FALSE) {
  #
  # Loads a GRangesList object with the exons grouped by gene and the
  #   names of the genes if it exists, otherwise it creates it
  #
  # Requires:
  #   data.table
  #   here
  #   GenomicFeatures
  #   biomaRt
  #   magrittr
  #   GenomicRanges
  #
  # Load command:
  # sapply(c("data.table", "here", "GenomicFeatures", "biomaRt", "magrittr",
  #       "GenomicRanges"), library, character.only = T)
  #
  # Args:
  #   TxDb:                     The database TxDb object
  #   organism:                 The name of the organism of interest e.g. "Homo sapies" or "Arabidopsis thaliana"
  #   genome.version:           The genome version to use
  #   genes.annotation.version: The gene annotation version to use
  #   protein.coding.only:      Default is FALSE Get the exons for the protein coding genes only
  #   overwrite:                Default is FALSE. Should I overwrite the file?
  #
  # Returns:
  #   A GRangesList object with the exons grouped by gene and the
  #   names of the genes
  #
  
  # Read the TxDb object for e.g. hg19 v19
  # Beware that chomosomes are stored as 
  # chr1, chr2, ..., chrX, chrY, chrMT
  
  # Prepare an empty exons granges object
  exons.granges <- NULL
  
  # Add the gencode prefix for human and mouse
  if(organism == "Homo sapiens" |
     organism == "Mus musculus") {
    
  }
  
  # Prepare the file path for the exons file
  exons.file <- file.path(here(),
                      "data-input",
                      "exon-coordinates",
                      paste0(genome.version, "-", genes.annotation.version, "-",
                             ifelse(protein.coding.only == TRUE,
                                    "exons-protein-coding-only.rds",
                                    "exons.rds")))
  
  # Prepare the exons bed file path
  exons.bed.file <- file.path(here(),
                          "data-input",
                          "BED-files",
                          paste0(genome.version, "-", genes.annotation.version, "-",
                                 ifelse(protein.coding.only == TRUE,
                                        "exons-protein-coding-only.bed",
                                        "exons.bed")))
  
  # Prepare the file path for the single exon genes file
  single.exons.file <- file.path(here(),
                             "data-input",
                             "single-exon-coordinates",
                             paste0(genome.version, "-", genes.annotation.version, "-",
                                    ifelse(protein.coding.only == TRUE,
                                           "single-exons-protein-coding-only.rds",
                                           "single-exons.rds")))
  
  # Prepare the file path for the single exon genes bed file
  single.exons.bed.file <- file.path(here(),
                                 "data-input",
                                 "BED-files",
                                 paste0(genome.version, "-", genes.annotation.version, "-",
                                        ifelse(protein.coding.only == TRUE,
                                               "single-exons-protein-coding-only.bed",
                                               "single-exons.bed")))
  
  # If the exons file does not exist, create it, otherwise read it
  if(file.exists(exons.file) == FALSE | overwrite == TRUE) {
    
    message("Get the exon coordinates...")
    
    # Now subset the genes to get only the common between genome annotation/ensembl
    exons.by.genes <- get.raw.exons.per.gene(TxDb,
                                             organism,
                                             genome.version,
                                             genes.annotation.version,
                                             protein.coding.only,
                                             overwrite = overwrite)
    
    message("Flank the coordinates on each side by 10 bases...")
    
    # Now for each gene we are gonna flank the exons
    # by 10 base pairs to each direction
    exons.by.genes@unlistData@ranges <- exons.by.genes@unlistData@ranges + rep(10,
                                                                               length(exons.by.genes@unlistData@ranges))
    
    #   =====     ====  =
    #     ====    ====    =
    #   Merge Exons
    #   ======    ====  = =
    #
    message("Merge the overlapping exons...")
    
    # I will merge overlaping exons
    reduced.exons.by.genes <- reduce(exons.by.genes)
    
    # I will merge overlaping exons
    exons.by.genes.merged.exons <- reduce(exons.by.genes)
    
    # Update the gene names
    exons.by.gene.names <- names(exons.by.genes.merged.exons)
    
    
    # Find single exon genes
    exons.per.gene.per.strand <- table(strand(exons.by.genes.merged.exons))
    
    # Sum the 3 columns of the table
    exons.per.gene <- rowSums(exons.per.gene.per.strand)
    
    # Find the genes with 1 exon
    single.exon.gene.names <- names(exons.per.gene[exons.per.gene == 1])
    
    # Get a GRangesList subset with the single exon genes
    single.exon.genes <- exons.by.genes.merged.exons[exons.by.gene.names %in% single.exon.gene.names]
    
    # Get the exons by genes subset without the single exon genes
    exons.by.genes <- exons.by.genes.merged.exons[!exons.by.gene.names %in% single.exon.gene.names]
    
    message("Save the exon coordinates as RDS...")
    
    # Save the exons in a file
    saveRDS(exons.by.genes,
            exons.file)
    
    # Save the single exon genes in a file
    saveRDS(single.exon.genes,
            single.exons.file)
    
    # And prepare them for export
    exons.granges <- exons.by.genes
    
    message("Save the exon coordinates as BED...")
    
    # Export exons as bed file
    export(exons.granges, exons.bed.file)
    
    # Export exons as bed file
    export(single.exon.genes, single.exons.bed.file)
    
    
    message("Exons file created!\n")
    message("Single exon genes file created")
  } else {
    message("Exons file found! Loading...\n")
    
    # Read the exons file
    exons.granges <- readRDS(exons.file)
  }
  
  return (exons.granges)
}

get.introns.per.gene <- function(TxDb, exons.by.genes,
                                 organism,
                                 genome.version,
                                 genes.annotation.version,
                                 protein.coding.only = FALSE,
                                 number.of.threads = 5,
                                 overwrite = FALSE) {
  #
  # Loads a GRangesList object with the introns grouped by gene and the
  #   names of the genes if it exists, otherwise it creates it
  #
  # Requires:
  #   data.table
  #   here
  #   parallel
  #   GenomicFeatures
  #   GenomicRanges
  #
  # Load command:
  # sapply(c("data.table", "here", "parallel", "GenomicFeatures",
  #       "GenomicRanges"), library, character.only = T)
  #
  # Args:
  #   TxDb:                     The database object for hg19 gencode v19
  #   exons.by.genes:           The exons GRangesList object for hg19 gencode v19
  #                             with exons grouped by gene, flanked by 10 bases from
  #                             each side and merged when they overlap
  #   organism:                 The name of the organism of interest e.g. "Homo sapies" or "Arabidopsis thaliana"
  #   genome.version:           The genome version to use e.g. gencode.v38
  #   genes.annotation.version: The gene annotation version to use
  #   protein.coding.only:      Default is FALSE Get the exons for the protein coding genes only
  #   number.of.threads:        Default is 5. How many threads should I use
  #   overwrite:                Default is FALSE. Should I overwrite the file?
  #
  #
  # Returns:
  #   A GRangesList object with the introns grouped by gene and the
  #   names of the genes
  #
  
  # Add the gencode prefix for human and mouse
  if(organism == "Homo sapiens" |
     organism == "Mus musculus") {
    
  }
  
  # Create the file path for the output file
  introns.file <- file.path(here(),
                            "data-input",
                            "intron-coordinates",
                            paste0(genome.version, "-", genes.annotation.version, "-",
                                   ifelse(protein.coding.only == TRUE,
                                          "introns-protein-coding-only.rds",
                                          "introns.rds")))
  
  # Prepare the introns bed file path
  introns.bed.file <- file.path(here(),
                            "data-input",
                            "BED-files",
                            paste0(genome.version, "-", genes.annotation.version, "-",
                                   ifelse(protein.coding.only == TRUE,
                                          "introns-protein-coding-only.bed",
                                          "introns.bed")))
  
  # Prepare the introns for export
  introns.granges <- NULL
  
  # If the exons file does not exist, create it, otherwise read it
  if(file.exists(introns.file) == FALSE | overwrite == TRUE) {
    
    message("Introns files not found. Let's create it...\n")
    
    # Then we will get the ranges for each gene from gencode
    gencode.gene.ranges <- genes(TxDb,
                                 single.strand.genes.only = FALSE)
    
    # And i am gonna flank them by 10 as well
    gencode.gene.ranges@unlistData@ranges <- gencode.gene.ranges@unlistData@ranges + rep(10, 
                                                                                         length(gencode.gene.ranges@unlistData@ranges))
    # Get gencode gene names
    gencode.gene.ranges.names <- names(gencode.gene.ranges)
    
    # Get the hg19v19 gene names
    exons.by.gene.names <- names(exons.by.genes)
    
    # I will subset the gencode gene ranges for the common genes for which i have exonic
    # information
    gencode.gene.ranges <- gencode.gene.ranges[gencode.gene.ranges.names %in% exons.by.gene.names]
    
    # Make an empty Grangeslist object 
    introns.by.genes <- GRangesList()
    
    # Make sure that both GRangesLists have the same order
    if( length(exons.by.genes) != length(gencode.gene.ranges) |
        all(names(exons.by.genes) == names(gencode.gene.ranges)) == FALSE){
      message("Invalid GRangesList objects...\n")
    }
    
    # Create the introns ranges list using multiple threads
    introns.by.genes <- mcmapply(create.introns.coordinates,
                                 exons.by.genes,
                                 gencode.gene.ranges,
                                 mc.cores = number.of.threads)
    
    # Make the list as a GRangesList object
    introns.by.genes <- GRangesList(introns.by.genes)
    
    # Save the introns file
    saveRDS(introns.by.genes,
            introns.file)
    
    # Prepare the introns for export
    introns.granges <- introns.by.genes
    
    # Export introns as bed file
    export(introns.granges, introns.bed.file)
    
    message("Introns files created!\n")
  } else {
    
    message("Introns file found! Loading...\n")
    
    # Read the existing introns file
    introns.granges <- readRDS(introns.file)
  }
  
  return (introns.granges)
}

create.introns.coordinates <- function(exons.by.gene, gencode.gene.range) {
  #
  # Creates the intronic coordinates of a gene by taking the exonic coordinates
  # of that gene and the total coordinates of that gene as stated in gencode
  # and returns the gaps between the introns
  #
  # Requires:
  #   data.table
  #   here
  #   parallel
  #   GenomicFeatures
  #   GenomicRanges
  #
  # Load command:
  # sapply(c("data.table", "here", "parallel", "GenomicFeatures",
  #       "GenomicRanges"), library, character.only = T)
  #
  # Args:
  #   exons.by.genes:     The exon granges for a specific gene
  #   gencode.gene.range: The coordinates of that gene according to gencode
  #   
  # Returns:
  #   The intronic ranges for a gene
  #
  
  # Prepare the chromosome factor for the granges object
  exons.by.gene.chromosome <- as.character(exons.by.gene@seqnames@values)
  exons.by.gene.chromosome.levels <- levels(exons.by.gene@seqnames)
  
  # Prepare the new granges for the granges object
  new.ranges <- gaps(exons.by.gene@ranges,
                     start = gencode.gene.range@ranges@start,
                     end = gencode.gene.range@ranges@start + gencode.gene.range@ranges@width-1)
  
  # Calculate the length of the new Rle factor
  exons.by.gene.chromosome.length <- length(new.ranges)
  
  # Make the new Rle factor for chromosome
  new.seqnames <- Rle(factor(c(exons.by.gene.chromosome),
                             exons.by.gene.chromosome.levels),
                      exons.by.gene.chromosome.length)
  
  # The seqinfor of the exons.by.gene
  exons.seqinfo <- seqinfo(exons.by.gene)
  
  # Create the seqinfo object
  new.seqinfo <- Seqinfo(seqnames = exons.seqinfo@seqnames,
                         seqlengths = exons.seqinfo@seqlengths,
                         isCircular = exons.seqinfo@is_circular,
                         genome = exons.seqinfo@genome)
  
  # Prepare as well the new Rle factor for the strand as well
  exons.by.gene.strand <- as.character(exons.by.gene@strand@values)
  exons.by.gene.strand.level <- levels(exons.by.gene@strand)
  
  # Make the new Rle factor for strand as well
  new.strand <- Rle(factor(c(exons.by.gene.strand),
                           exons.by.gene.strand.level),
                    exons.by.gene.chromosome.length)
  
  # Create the new seqlengths
  new.seqlengths <- seqlengths(exons.by.gene)
  
  # Finally make the new Granges object
  new.granges.item <- GRanges(seqnames = new.seqnames,
                              ranges = new.ranges,
                              strand = new.strand,
                              seqinfo = new.seqinfo,
                              seqlengths = new.seqlengths)
  
  # Finally set the seqlevels style
  seqlevelsStyle(new.granges.item) <- seqlevelsStyle(exons.by.gene)
  
  return(new.granges.item)
}

get.gene.info <- function(genome.version = "hg19",
                             protein.coding.only = FALSE,
                             overwrite = FALSE) {
  #
  # Get the gene lengths for a genome verion
  #
  # Requires:
  #   data.table
  #   magrittr
  #   here
  #   biomaRt
  # 
  # Load command:
  #  sapply(c("data.table", "magrittr", "biomaRt", "here"),
  #        library, character.only = T)
  #
  # Args:
  #   genome.version:       Default is "hg19". The genome version of choice hg19, hg38
  #                         mm9, mm10, TAIR10
  #   protein.coding.only:  Default is FALSE. Should I get the protein coding only?
  #   overwrite:            Default is FALSE. Should I overwrite the results?
  # 
  # Returns:
  #   A table with the results ensembl search with the gene lengths
  # 
  
  # Prepare the gene-lengths paths
  gene.lengths.folder.path <- file.path(here(),
                                        "data-input",
                                        "gene-lengths")
  
  # Prepare the gene lengths life name
  gene.lengths.path <- file.path(gene.lengths.folder.path,
                                 paste0(genome.version, "-gene-lengths",
                                        ifelse(protein.coding.only == TRUE,
                                               "-protein-coding",
                                               ""),
                                        ".csv"))
  
  if(file.exists(gene.lengths.path) == TRUE & overwrite == FALSE) {
    
    message("Gene lengths file for ", genome.version, " found! Loading...")
    
    ensembl.results <- fread(gene.lengths.path)
    
  } else {
    
    message("Gene lengths file for ", genome.version, 
            ifelse(protein.coding.only == TRUE,
                   " for protein coding only",
                   " for all types of transcripts"),
            " not found! Let's create it!")
    
    # Read the table with the hosts
    organisms.to.hosts <- fread(here("data-input",
                                     "organisms-to-hosts.csv"))
    
    # Save again to avoid mixing
    genome.version.input <- genome.version
    
    # Get the hostname
    hostname <- organisms.to.hosts[genome.assembly == genome.version.input, host]
    
    # Get the hostname
    mart.to.access <- organisms.to.hosts[genome.assembly == genome.version.input, mart]
    
    # Prepare to connect to Ensembl
    ensembl <- useMart(mart.to.access, host = hostname)
    
    # Get the dataset of interest
    dataset.of.interest <- listDatasets(ensembl) %>% 
      as.data.table %>%
      .[version == organisms.to.hosts[genome.assembly == genome.version.input,
                                      genome.version],
        dataset]
    
    # Use the dataset
    ensembl.mart <-  useDataset(dataset.of.interest, mart = ensembl)
    
    if(genome.version.input == "hg19"| 
       genome.version.input == "hg38") {
      
      ensembl.columns <- c("ensembl_gene_id_version",
                           "ensembl_transcript_id_version",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "transcript_start",
                           "transcript_end",
                           "transcript_length",
                           "transcript_biotype",
                           "entrezgene_accession",
                           "entrezgene_id",
                           "hgnc_symbol",
                           "external_gene_name",
                           "description")
      
      ensembl.filters <- c("transcript_biotype","chromosome_name")
      
      column.order <- c("ensembl_gene_id_version", "ensembl_transcript_id_version",
                        "chromosome_name", "start_position", "end_position", "gene_length",
                        "transcript_start", "transcript_end", "transcript_length",
                        "transcript_biotype", "entrezgene_accession", "entrezgene_id",
                        "hgnc_symbol", "external_gene_name", "description")
      
    } else if(genome.version.input == "mm10") {
      
      ensembl.columns <- c("ensembl_gene_id_version",
                           "ensembl_transcript_id_version",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "transcript_start",
                           "transcript_end",
                           "transcript_length",
                           "transcript_biotype",
                           "entrezgene_accession",
                           "entrezgene_id",
                           "mgi_symbol",
                           "external_gene_name",
                           "description")
      
      ensembl.filters <- c("transcript_biotype","chromosome_name")
      
      column.order <- c("ensembl_gene_id_version", "ensembl_transcript_id_version",
                        "chromosome_name", "start_position", "end_position", "gene_length",
                        "transcript_start", "transcript_end", "transcript_length",
                        "transcript_biotype", "entrezgene_accession", "entrezgene_id",
                        "mgi_symbol", "external_gene_name", "description")
      
    } else if(genome.version.input == "mm9") {
      
      ensembl.columns <- c("ensembl_gene_id",
                           "ensembl_transcript_id",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "transcript_start",
                           "transcript_end",
                           # "transcript_length",
                           "transcript_biotype",
                           "entrezgene",
                           "mgi_symbol",
                           "external_gene_id",
                           "description")
      
      ensembl.filters <- c("biotype","chromosome_name")
      
      column.order <- c("ensembl_gene_id", "ensembl_transcript_id",
                        "chromosome_name", "start_position", "end_position", "gene_length",
                        "transcript_start", "transcript_end", "transcript_length",
                        "transcript_biotype", "entrezgene",
                        "mgi_symbol", "external_gene_id", "description")
      
    } else if(genome.version.input == "TAIR10") {
      
      ensembl.columns <- c("ensembl_gene_id",
                           "ensembl_transcript_id",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "transcript_start",
                           "transcript_end",
                           "transcript_length",
                           "transcript_biotype",
                           "entrezgene_accession",
                           "entrezgene_id",
                           "tair_symbol",
                           "external_gene_name",
                           "description")
      
      ensembl.filters <- c("transcript_biotype","chromosome_name")
      
      column.order <- c("ensembl_gene_id", "ensembl_transcript_id",
                        "chromosome_name", "start_position", "end_position", "gene_length",
                        "transcript_start", "transcript_end", "transcript_length",
                        "transcript_biotype", "entrezgene_accession", "entrezgene_id",
                        "tair_symbol", "external_gene_name", "description")
    } else {
      stop("Invalid genome version: ", genome.version.input)
    }
    
    # Get the data from Ensembl    
    ensembl.results <- data.table(getBM(ensembl.columns,
                                        # filters = c("chromosome_name"),
                                        # values = list(c(as.character(1:22),
                                        #                 "X", "Y", "MT",
                                        #                 "Mt", "Pt")),
                                        mart = ensembl.mart))
    
    
    # Should I get the protein coding only?
    if(protein.coding.only == TRUE) {
      
      ensembl.results <- ensembl.results[get(ensembl.filters[1]) == "protein_coding" &
                                         get(ensembl.filters[2]) %in% c(1:22,
                                                                "X", "Y", "MT",
                                                                "Mt", "Pt")]
      
    } else {
      ensembl.results <- ensembl.results[chromosome_name %in% c(1:22,
                                                                "X", "Y", "MT",
                                                                "Mt", "Pt")]
      
    }
    # Update the mitochondrial chromosomes
    ensembl.results[chromosome_name %in% c("MT", "Mt"), chromosome_name := "M"]
    
    # And the chloloplasts
    ensembl.results[chromosome_name %in% c("Pt"), chromosome_name := "C"]
    
    
    if(genome.version.input == "mm9") {
      # Calculate the gene lengths
      ensembl.results[, transcript_length := transcript_end - transcript_start]
    }
    
    # Calculate the gene lengths
    ensembl.results[, gene_length := end_position - start_position]
    
    # Set the order of the data
    setcolorder(ensembl.results, column.order)
    
    # Fix the description column,
    ensembl.results[, description := gsub(",", ";", description)]
    
    # If the gene lengths dir does not exist, create it
    if(dir.exists(gene.lengths.folder.path) == FALSE) {
      dir.create(gene.lengths.folder.path)
    }
    
    message("Storing ensemble data table under: ", gene.lengths.path)
    
    # Write the file with the results
    fwrite(ensembl.results, gene.lengths.path, quote = F, sep = ",")
    
  }
  
  return (ensembl.results)
}

# ==== TODO make it more generic ==== 
get.sample.table.from.dataset.annotation <- function( annotation.table, 
                                              analysis.path,
                                              conditions.column,
                                              sample.column,
                                              bam.files.directory = "",
                                              samples.to.discard.pattern = "",
                                              conditions.of.interest = "all",
                                              number.of.samples = 5,
                                              samples.file.name = "annotation-samples.csv") {
  #
  # Reads the annotation table of a generic dataset containing information for all 
  # the experiments, and draws N random samples for the conditions.column of interest
  #
  # Args:
  #   annotation.table:       A data.table containing the metadata of Prokisch samples
  #   analysis.path           The analysis path string
  #   conditions.column:      The name of the column containing the condition annotation
  #   sample.column:      The name of the column containing the sample names
  #   conditions.of.interest: Default all tissues. A vector containing the tissues of interest
  #   number.of.samples:      Default is 5. The number of random samples to draw or "all" to get 
  #                           all samples from a tissue
  #   gtex.samples.file.name: Default is "gtex-samples.csv". The file name with the
  #                           gtex random samples
  #
  # Returns:
  #   A data.table containing N random samples of each tissue/subtissue
  #
  
  annotation.samples.file.path <- file.path(here(),
                                        "data-input",
                                        analysis.path,
                                        samples.file.name)
  
  if(file.exists(annotation.samples.file.path) == FALSE) {
    
    message(annotation.samples.file.path," not found, it will be generated...\n")
    
    # If there is no regex pattern to discard
    # patch the regex with ^$ (otherwise everything is true, see default behaviour of grep)
    if(samples.to.discard.pattern == "") {
      samples.to.discard.pattern <- "^$"
    }
    
    # Prepare the sample table for export
    samples.table <- annotation.table[get(conditions.column) %in% conditions.of.interest &
                                      !grepl(samples.to.discard.pattern, FIBROBLAST_ID) &
                                        RNA_ID != "", 
                                      c("RNA_ID",
                                         "FIBROBLAST_ID",
                                         "PATIENT_ID",
                                         "GROWTH_MEDIUM",
                                         "IS_RNA_SEQ_STRANDED",
                                         "GENDER")]
    # Construct the paths for the bam files
    rna.bam.folders <- file.path(bam.files.directory,
                           samples.table$RNA_ID,
                           "RNAout/paired-endout")
    
    # Check that all the paths exist
    if(all(dir.exists(rna.bam.folders)) == FALSE) {
      stop("Some bam files paths do not exist!\n")
    }
    
    bam.files <- sapply(rna.bam.folders, list.files, ".*\\.(sort\\.bam)$", full.names = TRUE)
    
    samples.table[, RNA_BAM_FILE := bam.files]
  
    message("Write the ", samples.file.name," file...\n")
    
    # Write the table with the SRA ids of interest
    fwrite(samples.table,
                annotation.samples.file.path,
                sep = ",")
    
  } else {
    
    message("Read the samples csv file.\n")
    
    # Read the gtex samples file
    samples.table <- fread(annotation.samples.file.path)
  }
  
  return (samples.table)
}


# ==== Obsolete functions ====

get.exons.per.gene.unmerged <- function(TxDb.hg19.gencode19) {
  #
  # Loads a GRangesList object with the exons grouped by gene but without flanking and merging the exons
  # and the names of the genes if it exists, otherwise it creates it
  #
  # Args:
  #   TxDb.hg19.gencode19: The database object for hg19 gencode v19
  #
  # Returns:
  #   A GRangesList object with the exons grouped by gene and the
  #   names of the genes
  #
  
  # Read the TxDb object for hg19 v19
  # Beware that chomosomes are stored as 
  # chr1, chr2, ..., chrX, chrY, chrMT
  
  # Prepare an empty exons granges object
  exons.granges.unmerged <- NULL
  
  # Prepare the file path for the exons file
  exons.unmerged.file <- file.path(here(),
                                   "data-input",
                                  "exons-unmerged.rds")
  
  # Prepare the exons bed file path
  exons.unmerged.bed.file <- file.path(here(),
                                   "data-input",
                                  "exons-unmerged.bed")
  
  # If the exons.unmerged file does not exist, create it, otherwise read it
  if(file.exists(exons.unmerged.file) == FALSE) {
    
    message("exons.unmerged file not found. Let's create it...\n")
    # Group the exons.unmerged by gene
    exons.unmerged.by.genes <- exonsBy(TxDb.hg19.gencode19, use.names= TRUE)
    
    # Prepare to connect to Ensembl for grch37
    ensembl <- useMart("ensembl", host = "grch37.ensembl.org")
    
    # Pick the human dataset
    ensembl.hg19 <-  useDataset("hsapiens_gene_ensembl", mart = ensembl)
    
    # Get the the information for all genes, all chromosomes as well the
    # type of transcript (protein coding or whatever)
    # and the description of the gene
    ensembl.genes <- data.table(getBM(attributes = c("ensembl_gene_id_version",
                                                     "chromosome_name",
                                                     "description",
                                                     "transcript_biotype",
                                                     "ensembl_transcript_id_version"),
                                      mart = ensembl.hg19))
    
    # Make a vector with all the standard chromosomes
    # In Ensembl chromosomes are stored as 1, 2, ..., X, Y, MT
    standard.chromosomes <- c( as.character(c(1:22)),
                               "X", 
                               "Y",
                               "MT")
    
    # Now subset the Ensembl data table by keeping only the protein coding
    # genes, and only those which belong to a standard chromosome
    ensembl.protein.coding.genes <- ensembl.genes[ transcript_biotype == "protein_coding" &
                                                     chromosome_name %in% standard.chromosomes, ]
    
    
    # Get the names of the protein coding genes for ensemble
    ensembl.protein.coding.gene.names <- ensembl.protein.coding.genes$ensembl_transcript_id
    
    # Get the hg19v19 gene names
    exons.unmerged.by.gene.names <- names(exons.unmerged.by.genes)
    
    # Now subset the gencode genes to get only the common between gencode/ensembl
    exons.unmerged.by.genes <- exons.unmerged.by.genes[exons.unmerged.by.gene.names %in% ensembl.protein.coding.gene.names]
    
    
    # Save the exons.unmerged in a file
    saveRDS(exons.unmerged.by.genes,
            exons.unmerged.file)
    
    
    # And prepare them for export
    exons.granges.unmerged <- exons.unmerged.by.genes
    
    message("Exons unmerged file created!\n")
  } else {
    message("Exons file found! Loading...\n")
    
    # Read the exons file
    exons.granges.unmerged <- readRDS(exons.unmerged.file)
  }
  
  return (exons.granges.unmerged)
}


