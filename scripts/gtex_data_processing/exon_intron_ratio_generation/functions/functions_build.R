# ==== Build functions ====

get.count.table.on.sample <- function(sample.name,
                                      analysis.name,
                                      features,
                                      bam.file,
                                      is.paired.end,
                                      count.fragments,
                                      sample.counts.output.dir) {
  
  #
  # Get the count table for one sample so that the pipeline is more modular
  #
  # Requires:
  #   data.table
  #   GenomicAlignments
  # 
  # Load command:
  #  sapply(c("data.table", "GenomicAlignments"),
  #        library, character.only = T)
  #
  # Args:
  #   sample.name:              The sample's name
  #   analysis.name:            The name for the folder of the analysis
  #   features:                 The features GRanges
  #   bam.file:                 The Bamfile of interest
  #   is.paired.end:            Default is TRUE. Are the experiments paired ended or not?
  #   count.fragments:          Default is FALSE. How to do the counting of the paire
  #   sample.counts.output.dir: The output dir for the sample RDS
  #
  # Returns:
  #   A count table with the number of reads that map to each feature (gene)
  #
  
  # Prepare the output folder
  sample.counts.output.path <- file.path(sample.counts.output.dir,
                                        paste0(sample.name, ".rds"))
  
  # If the sample is counted read it, otherwise count the bam file
  if(file.exists(sample.counts.output.path) == TRUE) {
    
    cat("Reading sample ", sample.name)
    
    sample.table <-  readRDS(sample.counts.output.path)
  } else {
    cat("Counting sample ", sample.name)
    
    retry(# IntersectionStrict demands from a  read to fall completely inside a feature
      {
        # Count the Bam File
        sample.table <- summarizeOverlaps(features = features,
                                          reads = bam.file,
                                          mode = "IntersectionStrict",
                                          singleEnd = !is.paired.end,
                                          ignore.strand = TRUE,   # Ignore the strand when counting
                                          inter.feature = FALSE,  # On FALSE counts mapping to multiple
                                          # features are counted once for each
                                          # feature
                                          fragments = count.fragments)
        },
      when = "Error",
    )
    
    # Save the sample's RDS
    saveRDS(sample.table, file = sample.counts.output.path)
  }
  
  return (list(sample.table))
}

get.counts.table <- function(features, analysis.name, table.type,
                             samples.to.condition, condition.column,
                             bam.file.list,
                             is.paired.end = TRUE, count.fragments = FALSE,
                             protein.coding.only = FALSE,
                             number.of.threads = 16, overwrite = FALSE) {
  #
  # Gets the features (exon or intron ranges for genes) and the bam files of the reads
  # and calculates the count table
  #
  # Requires:
  #   data.table
  #   GenomicAlignments
  # 
  # Load command:
  #  sapply(c("data.table", "GenomicAlignments"),
  #        library, character.only = T)
  #
  # Args:
  #   features:               The GRangesList object with the ranges of exons or introns
  #   analysis.name           The name for the folder of the analysis
  #   table.type:             Should be "exons", "introns" or "raw-exons"
  #   samples.to.condition:   A data.table with 4 columns mapping a bam file to a condition
  #   condition.column:       The column name with the conditions of interest
  #   bam.file.list:          The BamFilesList with the bam files
  #   is.paired.end:          Default is TRUE. Are the experiments paired ended or not?
  #   count.fragments:        Default is FALSE. How to do the counting of the paired ends.
  #   protein.coding.only:    Default is FALSE. Should i get only the coding genes
  #   number.of.threads:      Default 16. The number of processes to spawn
  #   overwrite:              Default is FALSE. Should I overwrite the results?
  #
  # Returns:
  #   A count table with the number of reads that map to each feature (gene)
  #
  
  # Create an empty count table 
  counts.table <- NULL
  
  # Prepare the RDS file name
  rds.file.name <- paste0(table.type,
                          "-counts-",
                          analysis.name,
                          ifelse(protein.coding.only == TRUE,
                                 "-protein-coding-only",
                                 ""),
                          ".rds")
  
  # Prepate the file name for the count table object
  counts.file.path <- here("data-output",
                           analysis.name,
                           "counts-rds",
                           rds.file.name)
  
  # Chech of the count table file exists, otherwise create it
  if(file.exists(counts.file.path) == FALSE | overwrite == TRUE) {
    # tryCatch({
      message("Count table for ", rds.file.name, " not found.Let's create it...")
      
      message(parallel::detectCores(), " cores were detected.")
      
      message("Counting ", table.type, " with ", number.of.threads, " threads...")
      
      sample.counts.output.dir <- file.path(here(),
                                            "data-output",
                                            analysis.name,
                                            "sample-counts-rds",
                                            table.type)
      
      # Create the analysis output dir if it does not exist
      if(dir.exists(sample.counts.output.dir) == FALSE) {
        
        message("Creating ", sample.counts.output.dir," dir...")
        
        dir.create(sample.counts.output.dir, recursive = TRUE)
      }
      
      # Try to calculate the samples
      retry(sample.counts <- pbmcmapply(get.count.table.on.sample,
                                      sample.name = names(bam.file.list),
                                      bam.file = bam.file.list,
                                      MoreArgs = list(is.paired.end = is.paired.end,
                                                      count.fragments = count.fragments,
                                                      features = features,
                                                      sample.counts.output.dir = sample.counts.output.dir,
                                                      analysis.name = analysis.name),
                                      mc.cores = number.of.threads,
                                      ignore.interactive = TRUE),
            when = "Error")
      
      # Prepare the multicore params
      # multicores.params <- MulticoreParam(workers = number.of.threads,
      #                                     progressbar = TRUE)
      
      counts.table <- Reduce(cbind, sample.counts)
      # IntersectionStrict demands from a  read to fall completely inside a feature
      # counts.table <- summarizeOverlaps(features = features,
      #                                   reads = bam.file.list,
      #                                   mode = "IntersectionStrict",
      #                                   singleEnd = !is.paired.end,
      #                                   ignore.strand = TRUE,   # Ignore the strand when counting
      #                                   inter.feature = FALSE,  # On FALSE counts mapping to multiple
      #                                                           # features are counted once for each
      #                                                           # feature
      #                                   fragments = count.fragments,
      #                                   BPPARAM = multicores.params)
      
      message(Sys.time(), " summarizeOverlaps for ", table.type, " done")
      
      # Create the path for the analysis output dir
      analysis.output.dir <- here("data-output",
                                  analysis.name)
      
      # Create the analysis output dir if it does not exist
      if(dir.exists(analysis.output.dir) == FALSE) {
          
        message("Creating ", analysis.output.dir," dir...")
        
        dir.create(analysis.output.dir, recursive = TRUE)
      }
      
      # Create an intermediate folder path
      intermediate.folder.path <- file.path(analysis.output.dir,
                                        "intermediate-tables")
      
      # Create an intermediate folder path
      counts.folder.path <- file.path(analysis.output.dir,
                                  "counts-rds")
      
      # Create the dir if it does not exist
      if(dir.exists(intermediate.folder.path) == FALSE) {
        dir.create(intermediate.folder.path, recursive = TRUE)
        dir.create(counts.folder.path, recursive = TRUE)
      }
      
      # Prepare the path for the intermediate results
      intermediate.file.path <- file.path(intermediate.folder.path,
                                      paste0("intermediate-", rds.file.name))
      
      # Save intermediate results
      saveRDS(counts.table, intermediate.file.path)
      
      # Store the old colData
      old.colData <- colData(counts.table)
  
      # Get the new colData bases on the file names
      new.colData <- generate.colData.frame(table.type,
                                            old.colData, 
                                            samples.to.condition,
                                              condition.column)
  
      # Update the new colData of the SummarizedExperiment
      colData(counts.table) <- new.colData
      
      # Save the count table to an RDS file
      saveRDS(counts.table, counts.file.path)
      
      message("Count table for ", rds.file.name, " created!")
    
      # },
      # error = function(e) { 
      #   stop(Sys.time()," Error in get.counts.table: ", conditionMessage(e), "\n")
      #   })
      
  } else {
    message("Count table for ", rds.file.name, " found! Loading...")
    
    # Load the count table
    counts.table <- readRDS(counts.file.path)
  }
  
  return (counts.table)
}

generate.colData.frame <- function(table.type, old.colData, samples.to.condition,
                                   condition.column) {
  #
  # Fills the colData dataframe giving each sample a type (exon/intron/plain-exon) and the 
  # condition on which belongs
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   table.type:           Should be "exons", "introns" or "raw-exons"
  #   old.colData:          The empty data.frame of colData containing only the
  #                         sample names as rownames
  #   samples.to.condition: The data.table with 2 columns, the RNA_ID and the condition
  #   condition.column:    The name of the condition column e.g. tissue  
  #
  # Returns:
  #   The updated colData data.frame filled the the origin of count (exon/intron)
  #   and the tissue
  #

  # Get the table type
  region.type <- table.type

  # Get the sample names
  sample.names <- rownames(old.colData)

  # Get the positions of the samples from the  samples.to.condition table
  # by the order they appear on the rownames
  condition.positions <- sapply(sample.names,
                                function(x, samples) grep(x, samples)[1],
                                samples.to.condition[, RNA_ID],
                                USE.NAMES = FALSE)

  # Now get the conditions on that order
  conditions <- samples.to.condition[condition.positions, get(condition.column)]

  # Get the number of samples
  number.of.samples <- dim(old.colData)[1]

  # Turn the tissues to factor and add the on the colData data.frame
  old.colData[[condition.column]] <- factor(conditions)

  # Add the region column with the exons/introns string
  old.colData[["region.type"]] <- factor(rep.int(region.type, number.of.samples))

  return (old.colData)
}

fix.SummarizedExperiment.condition.names <- function(count.table,
                                                      condition.column) {
  
  #
  # Fixes the condition names of a count table by replacing the "weird" characters
  # of the condition name ("(", "-", ")" etc.) with "_"
  # 
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   count.table:      The original data.table
  #   condition.column: The name of the condition column
  #
  # Returns:
  #   The fixed data.table with the altered levels for the tissues factors
  #
  
  # Fix tissue factors
  levels(colData(count.table)[[condition.column]]) <- gsub("(\ -\ )|(\ \\()|\ |-", 
                                                    "_",
                                                    levels(colData(count.table)[[condition.column]]))
  
  levels(colData(count.table)[[condition.column]]) <- gsub(")", 
                                                    "",
                                                    levels(colData(count.table)[[condition.column]]))
  
  return (count.table)
  
}

generate.tpms <- function(raw.exon.counts, analysis.name,
                          tpm.file.name = "TPM-table.csv",
                          genome.version,
                          genes.annotation.version,
                          protein.coding.only = FALSE, 
                          overwrite = FALSE) {
  #
  # Create a table of TPMs in case that no TPM table data exist by using the RNA-seq counts
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   raw.exon.counts:      The data.table containing the raw exon counts
  #   analysis.name:        The analysis name
  #   tpm.file.name:        Default is "TPM-table.csv".
  #   gene.lengths.file:    Default is "data-input/gene-lengths/gene.lengths.csv". The file
  #                         with the gene lengths of the major isoforms for each sample
  #   overwrite:            Default is FALSE. Should I overwrite the previous pseudoTPM table?
  #
  # Returns:
  #   The TPMs table path
  #
  
  # Prepare the TPM file path
  tpm.file.path <- file.path(here(),
                             "data-output",
                             analysis.name,
                             tpm.file.name)
  
  # If the TPM file does not exist or should be overwritten, create it
  if(file.exists(tpm.file.path) == FALSE | overwrite == TRUE) {
    
    message("TPM file not found. Let's create a TPM file...")
    
    # Get the raw exon counts table
    raw.exon.counts <- store.raw.counts(NULL, "raw-exons", analysis.name, return = TRUE)
    
    # Prepare the raw exons coordinates path
    raw.exons.granges.path <- file.path(here(),
                                        "data-input",
                                        "exon-coordinates",
                                        "raw-exons",
                                        paste0(genome.version, "-", genes.annotation.version, "-",
                                               ifelse(protein.coding.only == TRUE,
                                                      "raw-exons-protein-coding-only.rds",
                                                      "raw-exons.rds")))
    # Read the raw exons coordinates
    raw.exons.granges <- readRDS(raw.exons.granges.path)
    
    # Get the length of the exons for each gene
    gene.lengths <- data.table(gene.id = names(sum(width(raw.exons.granges))),
                               size.b = sum(width(raw.exons.granges)))
    
    # Get the gene lengths in kilobases
    gene.lengths[, size.kb := size.b/10^3]
    
    # And remove the obsolete size.b column
    gene.lengths[, size.b := NULL]
    
    
    # Melt the summed counts table
    raw.exon.counts.melted <- melt(raw.exon.counts,
                                   id.vars = "gene.id",
                                   variable.name = "sample.id",
                                   value.name = "counts")
    
    # Add the gene lengths
    raw.exon.counts.annotated <- merge(gene.lengths,
                                       raw.exon.counts.melted,
                                        by = "gene.id")
    
    # # Update the gene.id column
    # raw.combined.counts.melted.annotated[, gene.id := gene.id.version]
    # 
    # # Update the transcript id column
    # setnames(raw.combined.counts.melted.annotated, "transcript.id.version", "transcript.id")
    # 
    # # Remove the old column
    # raw.combined.counts.melted.annotated[, gene.id.version := NULL]
    # 
    # # Calculate the RPK for each gene on each sample
    raw.exon.counts.annotated[, RPK := counts/size.kb]
    
    # Now calculate the scaling factor for each sample
    raw.exon.counts.annotated[, scaling.factor := sum(RPK)/10^6,
                                         by = "sample.id"]
    
    # And finally the TPM of each gene on each sample
    raw.exon.counts.annotated[, TPM := RPK/scaling.factor]
    
    # Subset the table to keep only the needed columns
    raw.exon.counts.annotated <- raw.exon.counts.annotated[, .SD,
                                                           .SDcols = c("gene.id",
                                                                       "sample.id", 
                                                                       "TPM")]
    
    # Transform the table into a wide format
    tpm.per.sample <- dcast(raw.exon.counts.annotated,
                            gene.id ~ sample.id,
                            value.var = "TPM")
    
    # Fix the gene.id and transcript.id names
    setnames(tpm.per.sample,
             c("gene.id"),
             c("gene_id"))
    
    message("Save TPM file...")
    
    # And save it
    fwrite(tpm.per.sample,
           tpm.file.path,
           sep = ",",
           quote = F)
  } else {
    # Read the TPM table
    message("TPM file found! Loading...")
    
    tpm.per.sample <- fread(tpm.file.path)
  }
  
  return(tpm.file.path)
}

generate.pseudo.tpms <- function(raw.exon.counts, raw.introns.counts, analysis.name,
                                  tpm.file.name = "pseudoTPM-table.csv",
                                  genome.version = "hg19",
                                  overwrite = FALSE) {
  #
  # Create a table of pseudoTPMs in case that no TPM table data exist by summing the 
  # raw exons and raw introns tables
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   raw.exon.counts:      The data.table containing the raw exon counts
  #   raw.introns.counts:   The data.table containing the raw intron counts
  #   analysis.name:        The analysis name
  #   tpm.file.name:        Default is "pseudoTPM-table.csv".
  #   gene.lengths.file:    Default is "data-input/gene-lengths/gene.lengths.csv". The file
  #                         with the gene lengths of the major isoforms for each sample
  #   overwrite:            Default is FALSE. Should I overwrite the previous pseudoTPM table?
  #
  # Returns:
  #   The pseudoTPMs table path
  #
  
  # Prepare the TPM file path
  tpm.file.path <- here("data-output",
                          analysis.name,
                          tpm.file.name)
  
  # If the TPM file does not exist or should be overwritten, create it
  if(file.exists(tpm.file.path) == FALSE | overwrite == TRUE) {
    
    message("TPM file not found. Let's create a pseudoTPM file...")
    
    # Get the gene lengths of the protein coding genes
    gene.lengths <- get.gene.info(genome.version)
    
    # Keep only the gene id and the size in bases
    gene.lengths <- gene.lengths[, .SD, .SDcols = c("ensembl_gene_id_version",
                                                    "ensembl_transcript_id_version",
                                                    "gene_length")]
    
    # Rename the columns
    setnames(gene.lengths,
             c("ensembl_gene_id_version", "ensembl_transcript_id_version", "gene_length"),
             c("gene_id", "transcript_id","size.b"))
    
    # Get the longest transcript
    gene.lengths <- gene.lengths[gene.lengths[, .I[which.max(size.b)],
                                              by = "gene_id"]$V1]
    
    
    # Get the gene lengths in kilobases
    gene.lengths[, size.kb := size.b/10^3]
    
    # And remove the obsolete size.b column
    gene.lengths[, size.b := NULL]
    
    # Assume that the total reads of the prokisch for one gene
    # are the sum of exonic only and intronic only mapped reads
    raw.counts.summed <- raw.exon.counts[, -c(1)] + raw.introns.counts[, -c(1)]
    
    # Add the gene.id column
    raw.counts.summed <- cbind(raw.exon.counts[, 1], raw.counts.summed)
    
    # Melt the summed counts table
    raw.counts.summed.melted <- melt(raw.counts.summed, 
                                     id.vars = "gene.id",
                                     variable.name = "sample.id",
                                     value.name = "counts")
    
    # Set names of the gene id column to version
    setnames(raw.counts.summed.melted, "gene.id", "gene.id.version")
    
    # Add the gene.id column
    raw.counts.summed.melted[, gene.id := gsub("\\.[0-9]+", "", gene.id.version)]
    
    # Set names of the gene id column to version
    setnames(gene.lengths, c("gene_id", "transcript_id"),
             c("gene.id.version", "transcript.id.version"))
    
    # Add the gene.id column
    gene.lengths[, gene.id := gsub("\\.[0-9]+", "", gene.id.version)]
    
    
    # Add the gene lengths
    raw.combined.counts.melted.annotated <- merge(gene.lengths[, .(gene.id,
                                                                   transcript.id.version,
                                                                   size.kb)],
                                                  raw.counts.summed.melted,
                                                  by = "gene.id")
    
    # Update the gene.id column
    raw.combined.counts.melted.annotated[, gene.id := gene.id.version]
    
    # Update the transcript id column
    setnames(raw.combined.counts.melted.annotated, "transcript.id.version", "transcript.id")
    
    # Remove the old column
    raw.combined.counts.melted.annotated[, gene.id.version := NULL]
    
    # Calculate the RPK for each gene on each sample
    raw.combined.counts.melted.annotated[, RPK := counts/size.kb]
    
    # Now calculate the scaling factor for each sample
    raw.combined.counts.melted.annotated[, scaling.factor := sum(RPK)/10^6,
                               by = "sample.id"]
    
    # And finally the TPM of each gene on each sample
    raw.combined.counts.melted.annotated[, TPM := RPK/scaling.factor]
    
    # Subset the table to keep only the needed columns
    raw.combined.counts.melted.annotated <- raw.combined.counts.melted.annotated[, .SD,
                                                             .SDcols = c("gene.id",
                                                              "transcript.id",
                                                              "sample.id", 
                                                              "TPM")]
    
      # Transform the table into a wide format
    tpm.per.sample <- dcast(raw.combined.counts.melted.annotated,
                            gene.id + transcript.id ~ sample.id,
                            value.var = "TPM")
    
    # Fix the gene.id and transcript.id names
    setnames(tpm.per.sample,
             c("gene.id",
               "transcript.id"),
             c("gene_id",
               "transcript_id"))
    
    message("Save pseudoTPM file...")
    
    # And save it
    fwrite(tpm.per.sample,
           tpm.file.path,
           sep = ",",
           quote = F)
  } else {
    # Read the TPM table
    message("TPM file found! Loading...")
    
    # tpm.per.sample <- fread(tpm.file.path)
  }
  
  return(tpm.file.path)
}

log2.features <- function(features, size.factors = NULL, pseudocount = 1) {
  #
  # Gets the normalized features of interest (exons or introns) and transform them
  # by adding +1 pseudocound qnd logging
  #
  # Requires:
  #   data.table
  #
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   features:     The deseq2 object of interest (exons or introns)
  #   size.factors: Default is NULL. A vector with the size factor normalization values
  #   pseudocount:  Default is 1. Pseudocount to add.
  #
  # Returns:
  #   A log transformed data.table of the normalized features of interest
  #
  
  # Get the counts for exons or introns and wrap them in a data.table
  features.datatable <- as.data.table(counts(features)) 
  
  # Put the gene names as a column
  features.datatable[, gene.id := rownames(counts(features))]
  
  sample.columns <- colnames(features.datatable) 
  
  # # If the user provided a size.factor normalization vector,
  # # normalize table by size.factors
  # if(is.null(size.factors) == FALSE) {
  #   
  #   features.datatable <- do.size.factor.normalization(features.datatable, size.factors)
  # }
  
  # Find the sample columns
  sample.columns <- grep("gene.id", colnames(features.datatable), value = TRUE, invert = TRUE)
  
  # Reorder the columns
  setcolorder(features.datatable, c("gene.id", sample.columns))
  
  message("Add the pseudocount and log2 the counts...")
  
  # Add the pseudocount and log2 transform the table
  features.datatable.logged <- features.datatable[, (sample.columns) := log2(.SD + pseudocount), .SDcols = sample.columns]
  
  return (features.datatable.logged)
}

store.raw.counts <- function(features, feature.type, analysis.name,
                             protein.coding.only = FALSE,
                             return = FALSE,
                             overwrite = FALSE) {
  
  #
  # Exports the raw counts of exons/introns
  #
  # Requires:
  #   data.table
  #   SummarizedExperiment
  # 
  # Load command:
  #  sapply(c("data.table", "SummarizedExperiment"),
  #        library, character.only = T)
  #
  # Args:
  #   features:           The Summarized Experiment object for exons or introns
  #   feature.type:       The type of raw counts. Should be "exons", "introns" or "raw-exons"
  #   analysis.name:      The name of the analysis
  #   protein.coding.only:Default is FALSE. Should i get only the coding genes
  #   return:             Default is FALSE. If true it returns the raw.counts table of interest
  #                       (exon or intron)
  #   overwrite:          Default FALSE. Should I overwrite the raw counts file?
  #
  # Returns:
  #   Exports or returns the raw counts 
  #
  
  # If the feature type is wrong kill the execution
  if(feature.type %in% c("exons", "introns", "raw-exons") == FALSE) {
    stop("store.raw.counts: wrong feature names\n")  
  }
  
  # Prepare the file name for row counts
  raw.counts.file <- paste0("raw-counts-", feature.type, 
                            ifelse(protein.coding.only == TRUE,
                                   "-protein-coding-only",
                                   ""), ".csv")
  
  # Prepare the folder path
  raw.counts.folder <- here("data-output",
                            analysis.name,
                            "counts-csv")
  
  # Prepare the path for row counts
  raw.counts.path <- here(raw.counts.folder,
                           raw.counts.file)
  
  # If the raw counts csv does not exist lets create it
  if(file.exists(raw.counts.path) == FALSE | overwrite == TRUE) {
    
    # Create the analysis output dir if it does not exist
    if(dir.exists(raw.counts.folder) == FALSE) {
      message("Creating ", raw.counts.folder," dir...")
      dir.create(raw.counts.folder, recursive = TRUE)
    }
    
    # Prepare the raw counts for export
    raw.counts <- data.table(assay(features))
    
    # Add the gene.id column 
    raw.counts <- cbind(gene.id = rownames(rowData(features)),
                        raw.counts)
    
    message("Save raw ", feature.type, " counts in: ", raw.counts.path)
    
    # Export the raw counts
    fwrite(raw.counts, raw.counts.path, 
                quote = FALSE, sep = ",")
    
    return (NULL)
  } else if(return == TRUE) {
    
    message("Reading raw counts for: ", feature.type)
    # Load the raw counts table 
    raw.counts <- fread(raw.counts.path)
     
    return (raw.counts)
  } else {
    #message(raw.counts.file, " already exists.")
    return (NULL)
  }
}


combine.exons.introns <- function(exons.deseq2, introns.deseq2, analysis.name,
                                  annotation.table.path,
                                  mask.unexpressed.genes = TRUE,
                                  center.genes = TRUE,
                                  mask.by = "TPM",
                                  file.name = "EI-ratios.csv",
                                  file.name.masked = "EI-ratios-masked.csv",
                                  file.name.centered = "EI-ratios-centered.csv",
                                  overwrite = FALSE,
                                  tpm.table.path = "",
                                  condition.column = "subtissue",
                                  number.of.threads = 16,
                                  ...) {
  
  # Combine the exons and introns normalized counts by substraction introns
  # from exons and normalizing by the average difference
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   exons.deseq2:         The exons DeSeq2 object
  #   introns.deseq2:       The introns DeSeq2 object
  #   analysis.name:        The analysis name
  #   annotation.table.path:The path to the annotation of the experiment
  #   mask.by:              Default is "TPM". Method to do the masking. If method is "TPM"
  #                         it masks the genes on a sample level using the tpms, otherwise
  #                         masks them using the condition using the "condition".
  #   file.name:            Default is "EI-ratios.csv". The file name of the E/I ratios file.
  #   file.name.masked:     Default is "EI-ratios-masked.csv". The file name of the masked
  #                         E/I ratios file.
  #   file.name.centered:   Default is "EI-ratios-centered.csv". The file name of the 
  #                         centered E/I ratios file.
  #   tpm.table.path:       Defaults is "". The path to the tpm table. If not provided,
  #                         I will create a pseudoTPM table
  #   overwrite:            Default is FALSE. Should I overwrite the table?              
  #   number.of.threads:    Default 16. The number of processes to spawn
  #
  # Returns:
  #   The combinded data.table of exons/introns
  #
  
  # Prepare the file names for centered and non centered tables
  ei.ratios.path <- file.path(here(),
                              "data-output",
                              analysis.name,
                              file.name)
  
  ei.ratios.masked.path <- file.path(here(),
                                     "data-output",
                                     analysis.name,
                                     file.name.masked)
  
  
  ei.ratios.centered.path <- file.path(here(),
                                       "data-output",
                                       analysis.name,
                                       file.name.centered)
  
  # Does the "EI-ratios-centered.csv" exist if not create it otherwise load it
  if(file.exists(ei.ratios.centered.path) == TRUE & overwrite == FALSE) {
    
    message(file.name.centered, " found! Loading...")
    
    combined.data.table <- fread(ei.ratios.centered.path)
    
  } else {
    
    message(file.name.centered, " not found. Let's create it!")
    
    if(file.exists(ei.ratios.path) == TRUE & overwrite == FALSE) {
      
      message(file.name, " found! Loading...")
      
      combined.data.table <- fread(ei.ratios.path)
    } else {
      
      message("Generate EI ratios file (not centered)...")
      
      message("Get logged table of exons...")
      
      # Log2 transform the exons
      exons.logged <- log2.features(exons.deseq2)
      
      message("Get logged table of introns...")
      
      # Log2 transform the introns
      introns.logged <- log2.features(introns.deseq2)
      
      # Get the sample columns
      sample.columns <- grep("gene.id", colnames(exons.logged), value = TRUE, invert = TRUE)
      
      message("Get log2(Exon/Intron) table...")
      
      #### Get log2 Exon Intron ratio ####
      
      # Substract log2(intros) counts from log2(exons) to get the ratio
      # log2(exons/introns) = log2(exons) - log2(introns)
      combined.data.table <- exons.logged[, .SD, .SDcols = sample.columns] - introns.logged[, .SD, .SDcols = sample.columns] 
      
      # Add the gene column
      combined.data.table[, gene.id := exons.logged[, gene.id]]
      
      # Reset the column order
      setcolorder(combined.data.table, c("gene.id", sample.columns))
      
      # https://www.nature.com/articles/s41467-017-00867-z#Sec2
      # Inference of RNA decay rate from transcriptional profiling highlights 
      # the regulatory programs of Alzheimer’s disease 
      message("Remove the bias term")
      
      # ==== Remove bias ====
      combined.data.table <-  remove.ei.bias(combined.data.table,
                                             introns.logged,
                                             annotation.table.path,
                                             condition.column)
      
      
      message("Store EI ratios file (not centered)...")
      
      # Remove redundant objects
      rm(exons.logged, introns.logged)
      
      # Export the log2(EI) ratios table
      fwrite(combined.data.table,
             ei.ratios.path,
             sep = ",",
             quote = FALSE)
    }
    
    # Should I mask the unexpressed genes?
    if(mask.unexpressed.genes == TRUE) {
      
      # Does the masked file exists already?
      if(file.exists(ei.ratios.masked.path) & overwrite == FALSE) {
        
        message(file.name.masked, "found! Loading...")
        
        combined.data.table <- fread(ei.ratios.masked.path)
        
      } else {
        
        # Should i mask by a TPM threshold?
        if(mask.by == "TPM") {
          
          # Get the gene.ids
          gene.ids <- combined.data.table[, gene.id]
          
          # Prepare the annotation columns
          annotation.table.columns <- c("RNA_ID", condition.column, "Sample_Name")
          
          # Get the sample columns
          sample.columns <- grep("gene.id", colnames(combined.data.table), value = TRUE, invert = TRUE)
          
          # Did the user provided a TPM table path?
          # If not we will attempt to make a pseudoTPM table
          # if(tpm.table.path == "") {

            # message("No TPM table path provided. Will attempt to make a TPM table from scratch...")
            #
            # message("Loading raw exon counts...")
            #
            # # Get the raw exon counts table
            # raw.exon.counts <- store.raw.counts(NULL, "raw-exons", analysis.name, return = TRUE)


            # # Get the raw exon counts table
            # raw.exon.counts <- store.raw.counts(NULL, "exons", analysis.name, return = TRUE)
            #
            # # Get the raw intron counts table
            # raw.introns.counts <- store.raw.counts(NULL, "introns", analysis.name, return = TRUE)
            #
            # # Generate the pseudo TPM table and get the path
            # tpm.table.path <- generate.pseudo.tpms(raw.exon.counts,
            #                                         raw.introns.counts,
            #                                         analysis.name,
            #                                         overwrite = overwrite,
            #                                        genome.version = genome.version)

          # }
          
          message("Find genes detected by pipeline...")
          
          # Get raw exon counts
          raw.exon.counts <- store.raw.counts(NULL, "exons", analysis.name, return = TRUE)
          
          message("Mark genes detected in exonic regions...")
          
          # Melt the raw exon counts
          raw.exon.counts.melted <- melt(raw.exon.counts,
                                         id.vars = "gene.id",
                                         variable.name = "sample.id",
                                         value.name = "exon.counts")
          # Remove redundant object
          rm(raw.exon.counts)
          
          # Initialize column for genes detected in exons
          raw.exon.counts.melted[, is.detected.in.exons  := TRUE]
          
          # Mark genes not detected in exons
          raw.exon.counts.melted[exon.counts == 0, is.detected.in.exons  := FALSE]
          
          # Remove redundant column
          raw.exon.counts.melted[, exon.counts := NULL]
          
          # Get raw intron counts
          raw.intron.counts <- store.raw.counts(NULL, "introns", analysis.name, return = TRUE)
          
          message("Mark genes detected in intronic regions...")
          
          # Melt the raw introns counts
          raw.intron.counts.melted <- melt(raw.intron.counts,
                                           id.vars = "gene.id",
                                           variable.name = "sample.id",
                                           value.name = "intron.counts")
          
          # Remove redundant object
          rm(raw.intron.counts)
          
          # Initialize column for genes detected in introns
          raw.intron.counts.melted[, is.detected.in.introns  := TRUE]
          
          # Mark genes not detected in introns
          raw.intron.counts.melted[intron.counts == 0, is.detected.in.introns  := FALSE]
          
          # Remove redundant column
          raw.intron.counts.melted[, intron.counts := NULL]
          
          message("Merge genes detected in exons and introns by the pipeline...")
          
          # Merge tables for detected genes by the pipeline
          detected.in.pipeline <- merge(raw.exon.counts.melted,
                                        raw.intron.counts.melted,
                                        by = c("gene.id", "sample.id"))
          
          # Remove redundant objects
          rm(raw.exon.counts.melted, raw.intron.counts.melted)
          
          # Initialize column for genes detected
          detected.in.pipeline[, gene.is.detected := TRUE]
          
          # Mark undetected genes
          detected.in.pipeline <- detected.in.pipeline[is.detected.in.exons == FALSE &
                                                         is.detected.in.introns == FALSE, 
                                                       gene.is.detected := FALSE ]
          
          # Remove the redundant columns
          detected.in.pipeline[, is.detected.in.exons := NULL]
          detected.in.pipeline[, is.detected.in.introns := NULL]
          
          message("Load expressed genes per sample table...")
          
          # ==== Get expressed genes by sample table ====
          # Get the expressed genes table
         expressed.genes <- expressed.genes.per.condition(tpm.table.path = tpm.table.path,
                                                           annotation.table.path = annotation.table.path,
                                                           gene.ids = gene.ids,
                                                           sample.columns = sample.columns,
                                                           analysis.name = analysis.name,
                                                           annotation.table.columns = annotation.table.columns,
                                                           condition.column = condition.column,
                                                           number.of.threads = number.of.threads, 
                                                           overwrite.expressed.genes.file = overwrite,
                                                           ...)
          
          # Filter our samples that are not part of the analysis
          expressed.genes <- expressed.genes[, .SD, .SDcols = c("gene.id",
                                                                sample.columns)]
          
          message("Find non-expressed genes...")
          
          # Melt expressed genes table
          expressed.genes.melted <- melt(expressed.genes,
                                         id.vars = "gene.id",
                                         variable.name = "sample.id",
                                         value.name = "is.expressed")
          
          # Remove redundant objects
          rm(expressed.genes)
             
          # Merge the two masking tables tables
          masking.table <- merge(expressed.genes.melted,
                                 detected.in.pipeline,
                                 by = c("gene.id", "sample.id"))
          
          # Remove redundant objects
          rm(expressed.genes.melted, detected.in.pipeline)
          
          # Get genes/samples to keep
          masking.table[, to.keep := is.expressed & gene.is.detected]
          
          # Remove the redundant columns
          masking.table[, is.expressed := NULL]
          masking.table[, gene.is.detected := NULL]

          # Melt ei.ratios
          combined.data.table.melted <- melt(combined.data.table,
                                             id.vars = "gene.id",
                                             variable.name = "sample.id",
                                             value.name = "ei.ratio")
          
          # Remove redundant objects
          rm(combined.data.table)
          
          message("Merge ei ratios with the expressed genes table...")
          
          # Merge the two melted tables
          combined.data.table.melted <- merge(combined.data.table.melted,
                                               masking.table,
                                               by = c("gene.id", "sample.id"))
          
          # Remove redundant objects
          rm(masking.table)
          
          message("Keep only the expressed genes and those counted by the pipeline...")
          
          # Keep only genes expressed or caught by the pipeline
          combined.data.table.melted <- combined.data.table.melted[to.keep == TRUE]
          
          message("Finalize masked log2(E/I)...")
          
          #### Mask not expressed genes ####
          # Transpose the table and fill with NAs
          combined.data.table <- dcast(combined.data.table.melted,
                                       gene.id ~ sample.id,
                                       value.var = "ei.ratio",
                                       fill = NA_real_)
          
          # Remove redundant objects
          rm(combined.data.table.melted)
          
        } else {
          message("Mask base on condition...")
        }
        message("Store masked EI ratios...")
        
        # Export the log2(EI) ratios masked table
        fwrite(combined.data.table, ei.ratios.masked.path, sep = ",",
               quote = FALSE)
      }
    }
      
    # Should I center the genes?
    if(center.genes == TRUE) {
      message("Get means of genes across samples...")
      
      # Melt the table
      combined.data.table.melted <- melt(combined.data.table,
                                         id.vars = "gene.id",
                                         variable.name = "sample.id",
                                         value.name = "ei.ratio")
      
      # Get the gene mean
      combined.data.table.melted[, gene.average := mean(ei.ratio, na.rm = TRUE), by = "gene.id"]
      
      #### Center genes ####
      
      message("Center EI ratios...")
      
      # Remove the mean from the ei.ratio
      combined.data.table.melted[, ei.ratio := ei.ratio - gene.average]
      
      # Turn the table into wide format
      combined.data.table <- dcast(combined.data.table.melted,
                                   gene.id + gene.average ~ sample.id,
                                   value.var = "ei.ratio")
      
      # Reorder the columns
      setcolorder(combined.data.table, c("gene.id", "gene.average", sample.columns))
      
      message("Generate EI ratios file (centered)...")
      
      # Export the log2(EI) ratios centered table
      fwrite(combined.data.table,
             ei.ratios.centered.path, sep = ",",
             quote = FALSE)
      
      message(file.name.centered, " was created successfully!")
    } else {
      message("No centering of the table")
    }
    
  }  
  return (combined.data.table)
}

get.samples.annotation <- function(features.deseq2,
                                   analysis.name,
                                   annotation.table.path,
                                   condition.column,
                                   overwrite = FALSE) {
  
  #
  # Makes an annotation dataframe for the samples of the initial analysis
  # with the columns sample.id, condition and number of samples
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   features.colData:     The data.frame from the exons/introns DeSeq2 object
  #   analysis.name:        The string of the analysis.name (we use it as folder name)
  #   annotation.table.path:The path to the initial annotation of the samples
  #   condition.column:     The condition columm name e.g. subtissue for gtex or whatever
  #   file.name:            Default is "condition-annotation.csv".
  #                         The name of the exported file
  #   overwrite:            Default FALSE. Should I overwrite the raw counts file?
  #
  # Returns:
  #   A data.table with the annotation of the samples
  #
  
  # If no condition column is provided just exit
  if(is.null(condition.column) == TRUE) {
    return ()
  }
  
  message("======== Generate sample annotations ===========")
  
  # Prepare the file path
  condition.annotation.file.path <- here("data-output",
                                         analysis.name,
                                         paste0(condition.column,
                                                "-annotation.csv"))
    
    
  # Check if the condition annotation file already does not exist
  # If it doesn't create it otherwise load it
  if(file.exists(condition.annotation.file.path) == FALSE | overwrite == TRUE) {
    
    message("Conditions annotation file for ", condition.column," not found. Let's create it!")
    
    # Make a data.table from the data.frame of the deseq2 colData
    condition.annotation <- data.table(sample.id = rownames(features.deseq2))
    
    # Read the original annotation file from GTEx
    annotation.table <- fread(annotation.table.path)
    
    # Set the optional columns
    optional.columns <- c("Sample_Name", "sex", "RIN", "INDIVIDUAL_ID")
    
    # Check if some optional annotation columns exist
    if(all(c("RIN", "sex", "Sample_Name", "INDIVIDUAL_ID") %in% colnames(annotation.table)) == FALSE) {
      
      # If not all of the optional columns exists, get the ones that exist
      optional.columns <- grep(paste("(", optional.columns, ")", sep = "", collapse = "|"),
                               colnames(annotation.table),
                               value = TRUE)
    }
    
    # For the samples of the analysis, get the RNA_ID the RIN and the tissue columns
    annotation.table <- annotation.table[RNA_ID %in% condition.annotation[, sample.id],
                                         .SD,
                                         .SDcols = c("RNA_ID", optional.columns, condition.column)]
    
    # And now merge the 2 data.tables by sample id
    condition.annotation <- merge(condition.annotation,
                                   annotation.table,
                                   by.x = "sample.id", by.y = "RNA_ID")
    
    # And reorder the columns of the data.table
    setcolorder(condition.annotation, c("sample.id", optional.columns, sort(condition.column)))
    
    # Add the number of samples per condition
    condition.annotation[, eval(paste0("samples.per.", condition.column)) := .N,
                          by = condition.column]
    
    # Fix the column names
    setnames(condition.annotation, colnames(condition.annotation),
             tolower(gsub("[_-]", "\\.", colnames(condition.annotation))))
    
    message("Saving conditions annotation file for ", condition.column, " under: ",
            condition.annotation.file.path)
    
    fwrite(condition.annotation,
           condition.annotation.file.path,
           sep = ",",
           quote = FALSE)
    
  } else {
    
    message("Condition annotation file for ", condition.column," found! Loading...")
    
    # Read the tissue annotation file
    condition.annotation <- fread(condition.annotation.file.path)
    
  }
  return (list(condition.annotation))
}


average.ei.per.condition <- function(ei.ratios.per.condition,
                                     valid.percentage = 0.2) {
  # Averages the centered ei.ratios in one condition if valid values 
  # surpass a threshold, otherwise sets that gene on that condition to NA
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   ei.ratios.per.condition:A vector the centered ei.ratios
  #   valid.percentage:       Default is 0.2. The lower percentage threshold for
  #                           valid values for a gene to be expressed in one condition
  #
  # Returns:
  #   The mean of the vector or NA
  #
  
  # Check that the ei.ratios.per.condition is not null
  if(is.null(ei.ratios.per.condition) == TRUE) {
    return (NA_real_)
  }
  
  # Get the length of the valid ei ratios values in one condition
  valid.ei.ratios.per.condition.length <- length(ei.ratios.per.condition[!is.na(ei.ratios.per.condition)])
  
  # And the total length of all ei ratios values in that condition
  ei.ratios.per.condition.length <- length(ei.ratios.per.condition)
  
  # Calculate the ratio on valid values
  percentage.of.valid <- valid.ei.ratios.per.condition.length / ei.ratios.per.condition.length
  
  # If the valid values are equal or more that the valid.percentage threshold
  if(!is.na(percentage.of.valid) & percentage.of.valid >= valid.percentage) {
    
    # Average the centered ei ratios removing the NAs
    return (mean(ei.ratios.per.condition, na.rm = TRUE))
  } else {
    # Or add NA as real (NAs by default are logical)
    return (NA_real_)
  }
}

export.average.EI.by.condition <- function(EI.ratios.centered, condition.annotation, sample.columns,
                                           condition.column, analysis.name,
                                           file.name.suffix = "EI-ratios-centered-per",
                                           valid.values.percentage = 0.2,
                                           overwrite = FALSE) {
  
  #
  # Returns and export the averaged by condition data.table 
  # (rows are the genes/columns are the conditions)
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   EI.ratios.centered:   The data.table with the E/I ratios 
  #                         (columns are the sample names/rows are the gene names)
  #   condition.annotation: The data.table with the annotation of the samples
  #   sample.columns:       A vector with the column names of the samples
  #   condition.column:     The condition colummn name e.g. subtissue for gtex or whatever
  #   analysis.name:        The string of the analysis.name (we use it as folder name)
  #   file.name:            Default is "EI-ratios-centered-per".
  #                         The name of the exported file
  #   overwrite:            Default FALSE. Should I overwrite the raw counts file?
  #
  # Returns:
  #   The averaged E/I ratios data.table by condition
  #
  
  file.name <- paste(file.name.suffix,
                     paste0(tolower(condition.column),".csv"),
                     sep = "-")
  
  # Prepare the file path
  file.name.path <- here("data-output",
                         analysis.name,
                         file.name)
  
  # If the E/I ratios file exists just load it otherwise create it
  if(file.exists(file.name.path) == TRUE & overwrite == FALSE) {
    
    message("Loading previously created ",file.name, " ...")
    
    # Load the previously created file with gene names as factors and
    # set gene.id to be the key
    data.averaged <- fread(file.name.path)
    
    message(file.name, " loaded!")
    
  } else {
    
    message(file.name, " not found. Lets create it...")
    
    # Melt the EI ratios 
    data.melted <- melt(EI.ratios.centered[, .SD, .SDcols = c(sample.columns,
                                                                  "gene.id")],
                            id.vars = "gene.id",
                            variable.name = "sample.id",
                            value.name = "ei.ratio")
    
    message("Add condition annotation column...")
    
    # Merge the transposed data.table with the annotation column by the sample.id column
    merged.data <- merge(data.melted,
                         condition.annotation[, .SD, .SDcols = c("sample.id", condition.column)],
                         by = "sample.id")
    
    message("Average samples across conditions...")
    
    
    # Finalize the data into a long format
    data.averaged <- dcast(merged.data,
                           as.formula(paste("gene.id", condition.column,sep = "~")),
                           value.var = "ei.ratio",
                           fun.aggregate = average.ei.per.condition,
                           valid.percentage = valid.values.percentage)
    
    message("Exporting the averaged E/I ratios csv in: ", file.name.path ," ...")
    
    # Export the data.table
    fwrite(data.averaged, file.name.path, 
                quote = FALSE, sep = ",")
    
    message("Export Complete!")
  }
  
  
  return (list(data.averaged))
}

expressed.genes.per.condition <- function(tpm.table.path,
                                          annotation.table.path,
                                          gene.ids,
                                          sample.columns,
                                          analysis.name,
                                          annotation.table.columns = NULL,
                                          condition.column = "subtissue",
                                          tpm.threshold = 1,
                                          expressed.genes.file.prefix = "expressed-genes-by-sample",
                                          overwrite.expressed.genes.file = FALSE,
                                          number.of.threads = 16, ...) {
  #
  # Gets the expressed genes on the samples of each condition
  # 
  # Requires:
  #   data.table
  #   parallel
  # 
  # Load command:
  #  sapply(c("data.table", "parallel"),
  #        library, character.only = T)
  #
  # Args:
  #   tpm.table.path:               The TPM table
  #   annotation.table.path:        The annotation of the TPM table
  #   gene.ids:                     A vector with the genes I should keep
  #   sample.columns:               A vector with the samples I should keep
  #   analysis.name:                The name of the analysis
  #   annotation.table.columns:     Default is NULL. A vector with specific column 
  #                                 names of the tpm.annotation table in case the table
  #                                 is big (will load only those). Currently set to get
  #                                 the RNA_ID, Sample_Name and the condition.columns
  #                                 e.g. tissue, subtissue etc
  #   condition.column:             Default is "sutissue". The condition column from
  #                                 the sample annotation table
  #   tpm.threshold:                Default is 1. The TPM threshold on which genes
  #                                 to be called as expressed or not
  #   expressed.genes.file.prefix:  Default is gtex-expressed-genes-by-sample.csv". The file name of the expressed or not genes
  #   overwrite:                    Default is FALSE. Should I overwrite the 
  #                                 preexisting file?
  #   number.of.threads:            Default 16. The number of processes to spawn
  #
  # Returns:
  #   A table with one column with the gene ids and the rest of the columns are the condition columns containing the
  #   transcript id of the major isoform
  #
  
  # Prepare the path expressed genes on the pipeline
  path.to.expressed.genes <- here("data-output",
                                  analysis.name,
                                  paste0(expressed.genes.file.prefix,
                                         "-threshold-",
                                         tpm.threshold,
                                         ".csv"))
  
  # Prepare the path of all the expressed genes
  path.to.expressed.genes.all <- here("data-output",
                                       analysis.name,
                                       paste0(expressed.genes.file.prefix,
                                              "-threshold-",
                                              tpm.threshold,
                                              "-all-genes",
                                              ".csv"))
  
  # Does the isoforms file exists?
  if(file.exists(path.to.expressed.genes) == TRUE &
     overwrite.expressed.genes.file == FALSE) {
    
    message("Expressed genes per sample table file found! Loading...")
    
    # If yes load it
    expressed.genes.per.sample.pipeline <- fread(path.to.expressed.genes)
    
  } else {
    
    message("Expressed genes per sample table file not found...Let's create it!")
    
    # Should I select specific columns from the annotation table
    if(is.null(annotation.table.columns) == TRUE) {
      
      tpm.annotation.table <- fread(annotation.table.path)
    
      } else {
      tpm.annotation.table <- fread(annotation.table.path,
                                    select = annotation.table.columns)
    }
    
    # Get only the existing samples
    tpm.annotation.table <- tpm.annotation.table[RNA_ID %in% sample.columns]
    
    # If we want major isoforms per tissue or subtissue choose the appropriate column
    conditions <- tpm.annotation.table[, sort(unique(get(condition.column)))]
    
    # If the "Sample_Name" column does not exist in the annotation
    # copy the RNA_ID column
    if("Sample_Name" %in% colnames(tpm.annotation.table) == FALSE) {
      
      tpm.annotation.table[, Sample_Name := RNA_ID]
    }
    
    # The samples to select
    sample.names.to.select <- tpm.annotation.table[, Sample_Name]
    
    # Get a list with the sample names on each condition
    samples.per.condition <- tpm.annotation.table[, .SD,
                                                  .SDcols = c(condition.column,
                                                              "Sample_Name")] %>%
                              split(., by = condition.column,
                                    keep.by = FALSE)
    
    message("Get the expressed genes per condition")
    
    # Now iterate through all the conditions to get the major isoforms for each
    
    # Beware #
    # On purpose single thread as it explodes the memory for GTEx
    # Find the expressed genes for all conditions
    expressed.genes.list <- mapply(gene.is.expressed.in.samples,
                                     condition = as.list(names(samples.per.condition)),
                                     sample.names.to.select = samples.per.condition,
                                     MoreArgs = list(tpm.table.path = tpm.table.path,
                                                     tpm.annotation.table = tpm.annotation.table,
                                                     gene.ids = NULL,
                                                     condition.column = condition.column,
                                                     tpm.threshold = tpm.threshold))
    
    message("Merge the expressed genes per condition")
    
    # Combine all the list items
    expressed.genes.per.sample <- Reduce(merge, expressed.genes.list)
    
    message("And export them into the ", path.to.expressed.genes)
    
    # And write the expressed genes tables
    fwrite(expressed.genes.per.sample, 
           path.to.expressed.genes.all,
           sep = ",", quote = FALSE)
    
    
    # Subset for the genes of the pipeline
    expressed.genes.per.sample.pipeline <- expressed.genes.per.sample[gene.id %in% gene.ids]
    
    # And write the expressed genes tables
    fwrite(expressed.genes.per.sample.pipeline, 
           path.to.expressed.genes,
           sep = ",", quote = FALSE)
  }
  
  return (expressed.genes.per.sample.pipeline)
}

gene.is.expressed.in.samples <- function(condition.of.interest, 
                                         tpm.table.path, tpm.annotation.table,
                                         condition.column, gene.ids = NULL, 
                                         sample.names.to.select = NULL,
                                         tpm.threshold = 1) {
  #
  # Build a table with the expressed genes on each condition based on a threshold.
  # If the table contains transcripts as well, it picks the major transcript
  #
  # Requires:
  #   data.table
  # 
  # Load command:
  #  sapply(c("data.table"),
  #        library, character.only = T)
  #
  # Args:
  #   condition.of.interest:  The condition of interest
  #   tpm.table.path:         The tpm table path   
  #   tpm.annotation.table:   The table containing the annotation of the samples
  #   condition.column:       The column name of the column containing the conditions
  #   gene.ids:               Default is NULL. The genes of interest
  #   sample.names.to.select: Default is NULL. The sample names to target
  #   tpm.threshold:          Default is 1. The threshold on which one gene should be
  #                           considered expressed
  #
  # Returns:
  #   A data.table containing the samples in the columns and genes in the rows with TRUE or
  #   false values depending on the fact that the genes is expressed or not
  #
  
  # By default assume that the table has no transcripts inside
  table.has.transcripts <- FALSE
  
  message("Get TPM table for ", condition.of.interest)
  
  # Is the tpm table path on gct format?
  format.is.gct <- grepl("(\\.gct)", tpm.table.path)
  
  # Read the tpm table of the condition of interest
  if(is.null(sample.names.to.select) == TRUE) {
    
    # If the tpm table is a gct file, skip the first 2 rows
    if(format.is.gct == TRUE) {
      rows.to.skip <- 2
    } else {
      rows.to.skip <- 0 
    }
    
    # Read the tpm table
    tpm.table <- fread(tpm.table.path,
                       skip = rows.to.skip)  
  } else {
    if(format.is.gct == TRUE) {
      rows.to.skip <- 2
    } else
      rows.to.skip <- 0
      
    }
    
    tpm.table <- fread(tpm.table.path, select = c("gene_id",
                                                "transcript_id",
                                                sample.names.to.select[, Sample_Name]),
                       skip = rows.to.skip)
  
  # Should I get all the genes or only some of them?
  if(is.null(gene.ids) == FALSE) {
    
    # Get only the taget genes
    tpm.table <- tpm.table[gene_id %in% gene.ids]
  }
  
  # Get the sample names for the target condition
  samples.names <- tpm.annotation.table[get(condition.column) == condition.of.interest,
                                        Sample_Name]
  
  # Does the tpm.table has a transcript_id column?
  if("transcript_id" %in% colnames(tpm.table)) {
    table.has.transcripts <- TRUE
    
    # If yes the ensemble columns will be the transcript_ids and the gene_ids
    ensembl.id.columns <- c("transcript_id", "gene_id")
  } else {
    
    # Otherwise will be the gene.ids only
    ensembl.id.columns <- c("gene_id")
  }
  
  # Subset the TPM table for the target samples and genes
  tpm.table.subset <- tpm.table[,.SD, .SDcols = c(ensembl.id.columns, samples.names)]
  
  # Melt the table
  tpm.table.subset.melted <- melt(tpm.table.subset,
                                  id.vars = ensembl.id.columns,
                                  variable.name = "sample.name",
                                  value.name = "tpm")
  
  # If we have transcripts, get the major isoform
  if(table.has.transcripts == TRUE) {

    # Get the major trasncript for each gene by taking the max
    # If you have a tie on some gene, get the first
    tpm.table.subset.melted <- tpm.table.subset.melted[tpm.table.subset.melted[, .I[which.max(tpm)],
                                                                               by = c("sample.name", "gene_id")][, V1]]

  }
  
  # Initialize the is.expressed column
  tpm.table.subset.melted[, is.expressed := FALSE]
  
  # If the TPM of a gene in one sample is higher than the threshold
  # Set that gene to be expressed
  tpm.table.subset.melted[tpm >= tpm.threshold, is.expressed := TRUE]
  
  # Order the ensembl.id.columns
  ensembl.id.columns <- sort(ensembl.id.columns)
  
  # Reset the order of the columns
  setcolorder(tpm.table.subset.melted,
              c(ensembl.id.columns, "sample.name", "is.expressed", "tpm"))
  
  # Correct the ensembl.id.columns by replacing "_" with "."
  corrected.ensembl.id.columns <- gsub("_", "\\.", ensembl.id.columns)
  
  # Rename the columns
  setnames(tpm.table.subset.melted,
           ensembl.id.columns,
           corrected.ensembl.id.columns)
  
  # Merge the TPM table with the annotation
  tpm.table.subset.melted.annotated <- merge(tpm.table.subset.melted,
                                             tpm.annotation.table[Sample_Name %in% samples.names,
                                                                  .SD, .SDcols = c("Sample_Name",
                                                                                  "RNA_ID")],
                                             by.x = "sample.name",
                                             by.y = "Sample_Name")
  
  # Subset the TPM table to get only the relevant columns
  tpm.table.subset.melted.annotated <- tpm.table.subset.melted.annotated[, .SD,
                                                                         .SDcols = c(corrected.ensembl.id.columns,
                                                                                      "RNA_ID",
                                                                                      "is.expressed")]
  
  # And finally turn in into a wide format
  tpm.dcasted <- dcast(tpm.table.subset.melted.annotated, 
                       gene.id ~ RNA_ID,
                       value.var = "is.expressed")
  
  return (list(tpm.dcasted))
}

get.average.tpm.per.condition <- function(tpm.table.path, condition.column,
                                          condition.annotation, analysis.name,
                                        overwrite = FALSE) {
  
  #
  # Get the average TPM per transcript
  #
  # Requires:
  #   data.table
  #   magrittr
  #   here
  # 
  # Load command:
  #  sapply(c("data.table", "magrittr", "here"),
  #        library, character.only = T)
  #
  # Args:
  #   tpm.table.path:       The path to the TPM table file
  #   condition.column:     The name of the condition column
  #   condition.annotation: The condition annotation table 
  #   analysis.name:        The analysis name
  #   overwrite:            Default is FALSE. Should I overwrite the results?
  #
  # Returns:
  #   A table with the mean TPM per transcript
  # 
  
  if(is.null(condition.column) == TRUE) {
    return (list(NULL))  
  }
  
  # Prepare the tpm folder path
  average.tpm.folder.path <- here("data-output",
                                  analysis.name,
                                  "tpm-per-condition")
  
  # If the table folder does not exist, create it
  if(dir.exists(average.tpm.folder.path) == FALSE) {
    dir.create(average.tpm.folder.path, recursive = T)
  }
  
  # Prepare the table path
  average.tpm.path <- file.path(average.tpm.folder.path,
                            paste0("average-tpm-table-per-", condition.column, ".csv"))
  
  # Shoud I read or create the file
  if(file.exists(average.tpm.path) == TRUE & overwrite == FALSE) {
    
    message("Average TPM table for ", condition.column, " already exists! Loading...")
    
    # Read the table
    average.tpm.per.condition.wide <- fread(average.tpm.path)
  } else {
    message("Average TPM table for ", condition.column, " does not exist... Let's create it!")
    
    # Get the samples for each condition
    samples.tables <- condition.annotation %>%
      .[, .SD, .SDcols = c("sample.name", condition.column)] %>%
      split(., by = condition.column, keep.by = FALSE)
    
    # Get the average TPM on conditions
    average.tpm.on.conditions.list <- mapply(average.tpm.on.condition,
                                             condition.of.interest = names(samples.tables),
                                             samples.table = samples.tables,
                                             MoreArgs = list(tpm.table.path = tpm.table.path))
    
    # Merge the chunks
    average.tpm.per.condition <- rbindlist(average.tpm.on.conditions.list)
    
    # Turn the table into a wide format
    average.tpm.per.condition.wide <- dcast(average.tpm.per.condition,
                                            gene_id + transcript_id ~ condition,
                                            value.var = "mean.TPM")
    
    fwrite(average.tpm.per.condition.wide,
           average.tpm.path,
           quote = FALSE,
           sep = ",")
  }
  
  return (list(average.tpm.per.condition.wide))
}

get.tpm.per.condition <- function(tpm.table.path, analysis.name,
                                  condition.name,
                                  on.samples = FALSE,
                                  expressed.genes.by.sample.file = "expressed-genes-by-sample-threshold-1.csv") {
  
  #
  # Get the TPM per samples or conditions
  #
  # Requires:
  #   data.table
  #   here
  #
  # Load command:
  #   sapply(c("data.table", "here"),library, character.only = T)
  #
  # Args:
  #   tpm.table.path:                 The path to the TPM file
  #   analysis.name:                  The name of the analysis
  #   condition.name:                 The name of the condition
  #   on.samples:                     Default is FALSE. Should I keep the data on samples or 
  #                                   average the condition
  #   expressed.genes.by.sample.file: Default is "expressed-genes-by-sample-threshold-1.csv". The name of file
  #                                   for the expressed genes
  #
  # Returns:
  #   The table of the TPM per samples or averaged by conditions
  #
  
  # Prepare the folder path
  TPM.per.condition.folder.path <- file.path(here(),
                                             "data-output",
                                             analysis.name,
                                             "tpm-per-condition")
  
  # Create the folder if it does not exist
  if(dir.exists(TPM.per.condition.folder.path) == FALSE) {
    dir.create(TPM.per.condition.folder.path, recursive = TRUE)
  }
  
  # Prepare the file path
  TPM.per.condition.path <- file.path(TPM.per.condition.folder.path,
                                      paste0("tpm-table-per-",
                                             ifelse(on.samples == FALSE,
                                                    condition.name,
                                                    paste0("samples-by-", condition.name)),
                                             ".csv"))
  
  # If the table that I want exists, read it otherwise create it
  if(file.exists(TPM.per.condition.path) == TRUE) {
    
    TPM.per.condition <- fread(TPM.per.condition.path)
  } else {
    
    # Is the tpm table path on gct format?
    format.is.gct <- grepl("(\\.gct)", tpm.table.path)
    
    # If the tpm table is a gct file, skip the first 2 rows
    if(format.is.gct == TRUE) {
      rows.to.skip <- 2
    } else {
      rows.to.skip <- 0 
    }
    
    # Read the tpm table
    tpm.table <- fread(tpm.table.path,
                       skip = rows.to.skip)  
    
    # Prepare the annotation table file path
    condition.annotation.file.path <- file.path(here(),
                                                "data-output",
                                                analysis.name,
                                                paste0(condition.name, "-annotation.csv"))
    
    # Read and prepare the annotation table for the condition of interest
    annotation.table <- fread(condition.annotation.file.path) %>%
      .[, .SD, .SDcols = c("sample.id", condition.name, "sample.name")]
    
    # Prepare the expressed genes by sample file path
    expressed.genes.by.sample.file.path <- file.path(here(),
                                                     "data-output",
                                                     analysis.name,
                                                     expressed.genes.by.sample.file)
    
    # Read and prepare the expressed genes by sample
    expressed.genes.by.sample <- fread(expressed.genes.by.sample.file.path) %>%
      melt(.,
           id.vars = "gene.id",
           variable.name = "sample.id",
           value.name = "is.expressed")
    
    # Prepare the major isoforms per condition file path
    major.isoforms.per.condition.file.path <- file.path(here(),
                                                        "data-output",
                                                        analysis.name,
                                                        "major-isoforms",
                                                        paste0("major-isoforms-by-", condition.name, ".csv"))
    
    # Read and prepare the major isoforms per condition
    major.isoforms.per.condition <- fread(major.isoforms.per.condition.file.path) %>%
      melt(.,
           id.vars = "gene.id",
           variable.name = condition.name,
           value.name = "transcript.id")
    
    # Filter the TPM table for expressed genes and major isoforms
    tpm.table.filtered <- tpm.table %>%
      .[gene_id %in% expressed.genes.by.sample[, unique(gene.id)] &
          transcript_id %in% major.isoforms.per.condition[, unique(transcript.id)]] %>%
      melt(., 
           id.vars = c("gene_id", "transcript_id"),
           variable.name = "sample.name",
           value.name = "TPM")
    
    # Free up some space
    rm(tpm.table)
    invisible(gc())
    
    # Fix the names
    setnames(tpm.table.filtered, c("gene_id", "transcript_id"), c("gene.id", "transcript.id"))
    
    # Add the annotation
    tpm.table.filtered <- merge(tpm.table.filtered,
                                annotation.table,
                                by = "sample.name")
    
    # Add the isoforms to reduce the size table
    tpm.table.filtered <- merge(tpm.table.filtered,
                                major.isoforms.per.condition,
                                by = c("gene.id", condition.name, "transcript.id"))
    
    # Should I get the TPM per sample or condition
    if(on.samples == FALSE) {
      
      TPM.per.condition <- dcast(tpm.table.filtered,
                                 as.formula(paste0("gene.id ~ ", condition.name)),
                                 value.var = "TPM",
                                 fun = mean)
    } else {
      TPM.per.condition <- dcast(tpm.table.filtered,
                                 "gene.id ~ sample.id",
                                 value.var = "TPM")
    }
    
    # Save the file
    fwrite(TPM.per.condition,
           TPM.per.condition.path,
           quote = F, sep = ",")
  }
  
  return (TPM.per.condition)
}

remove.ei.bias <- function(combined.data.table,
                           introns.logged,
                           annotation.table.path,
                           condition.column) {
  
  # Remove the bias from the ei ratio
  # # https://www.nature.com/articles/s41467-017-00867-z
  # Inference of RNA decay rate from transcriptional profiling highlights 
  # the regulatory programs of Alzheimer’s disease 
  #
  # Requires:
  #   data.table
  #   fastglm
  # 
  # Load command:
  #  sapply(c("data.table", "fastglm"),
  #        library, character.only = T)
  #
  # Args:
  #   combined.data.table:  The table with the ei ratios       
  #   introns.logged:       The table with the logged intron counts 
  #   annotation.table.path:The path to the annotation of the experiment
  #   condition.column:     The name of the condition column
  #
  # Returns:
  #   The combinded data.table of exons/introns ratios with the removed bias
  #
  
  set.seed(24101992)
  
  # Melt the combined data table
  combined.data.table.melted <- melt(combined.data.table,
                                     id.vars = "gene.id",
                                     variable.name = "sample.id",
                                     value.name = "ei.ratio")
  
  # Read the annotation table and get the condition column and the sample.id
  condition.annotation <- fread(annotation.table.path, select = c("RNA_ID",
                                                                  condition.column))
  
  # Rename the column
  setnames(condition.annotation, "RNA_ID", "sample.id")
  
  # Add the annotation to the combined table
  combined.data.table.melted.annotated <- merge(combined.data.table.melted,
                                                condition.annotation,
                                                by = "sample.id")

  # Melt the intron counts  
  introns.logged.melted <- melt(introns.logged,
                                id.vars = "gene.id",
                                variable.name = "sample.id",
                                value.name = "introns.logged") 
  
  # Add the introns
  combined.data.table.melted.annotated <- merge(combined.data.table.melted.annotated,
                                                introns.logged.melted,
                                                by = c("sample.id", "gene.id"))
  
  # Fit a model for each gene for ei.ratio ~ introns.logged
  combined.data.table.melted.annotated[, fitted.introns.logged := fitted(fastglm(x = matrix(introns.logged,
                                                                                            ncol = 1),
                                                                                 y = ei.ratio)),
                                       by = c(condition.column, "gene.id")]
  
  # Remove the bias
  combined.data.table.melted.annotated[, ei.ratio.no.bias := ei.ratio - fitted.introns.logged]
  
  # Finalize the bias corrected table
  combined.data.table.no.bias <- dcast(combined.data.table.melted.annotated,
                                       gene.id ~ sample.id,
                                       value.var = "ei.ratio.no.bias")
  
  return (combined.data.table.no.bias)
  
}

# ==== Old ====
average.tpm.on.condition <- function(condition.of.interest, samples.table, tpm.table.path) {
  
  #
  # Get the average TPM per transcript on a condition
  #
  # Requires:
  #   data.table
  #   magrittr
  # 
  # Load command:
  #  sapply(c("data.table", "magrittr"),
  #        library, character.only = T)
  #
  # Args:
  #   condition.of.interest:  The condition of interest
  #   samples.table:          The table with the sample names on that condition
  #   tpm.table.path:         The path to the TPM table file
  #
  # Returns:
  #   A table with the mean TPM per transcript
  # 
  
  message("Getting the mean TPM on: ", condition.of.interest)
  
  # Read the tpm table for the samples of interest
  tpm.table <- fread(tpm.table.path,
                     select = c("gene_id", 
                                "transcript_id",
                                samples.table[, sample.name])) 
  
  # If transcript_id column does not exist, use duplicate the gene.id column
  if(("transcript_id" %in% colnames(tpm.table)) == FALSE) {
    tpm.table[, transcript_id := gene_id]
  }
  
  # Turn it in a wide format
  tpm.table.wide <- melt(tpm.table,
                         id.vars = c("gene_id", "transcript_id"),
                         variable.name = "sample.name",
                         value.name = "TPM")
  
  # Get the mean TPM per transcript
  tpm.table.wide[, mean.TPM := mean(TPM),
                 by = c("gene_id",
                        "transcript_id")] 
  
  # Get the unique values
  tpm.table.wide <- tpm.table.wide[, unique(.SD), .SDcols = c("gene_id",
                                                              "transcript_id",
                                                              "mean.TPM")]
  
  # Add the condition   
  tpm.table.wide[, condition := condition.of.interest]
  
  return (list(tpm.table.wide))
}

# ==== Dirty functions ====
# ==== Obsolete functions ====

do.size.factor.normalization <- function(features.datatable, size.factors) {
  
  #
  # Normalizes a count table by size factors
  #
  # Args:
  #   features.datatable: A count data table either with exons or introns
  #   size.factors:       A vector with the estimated size factors
  #
  # Returns:
  #   The normalized count table by size factors
  #
  message("Prepare to do size.factor normalization.")
  
  # Transpose the count table sample X genes
  features.datatable.transpose <- dcast(melt(features.datatable,
                                             id.vars = "gene.id",
                                             variable.name = "sample",
                                             value.name = "counts"),
                                        sample~gene.id,
                                        value.var = "counts")
  
  # Get the gene column names
  gene.columns <- grep("ENSG", colnames(features.datatable.transpose), value = TRUE)
  
  message("Do size.factor normalization.")
  
  # Divide each gene on all samples with the same factor
  features.datatable.transpose.normalized <- features.datatable.transpose[, mclapply(.SD, "/", size.factors),
                                                                          .SDcols = gene.columns]
  
  # Add the sample column
  features.datatable.transpose.normalized[, sample:= features.datatable.transpose$sample]
  
  # Reorder the table
  setcolorder(features.datatable.transpose.normalized, c("sample", gene.columns))
  
  # Transpose again the normalized table to gene.id X samples
  features.datatable <- dcast(melt(features.datatable.transpose.normalized,
                                   id.vars = "sample",
                                   variable.name = "gene.id",
                                   value.name = "normalized.counts"),
                              gene.id ~ sample,
                              value.var = "normalized.counts")
  
  return (features.datatable)
}

# get.gene.lengths <- function(features, region.type,
#                              dataset.name = "gtex") {
#   #
#   # Get the gene lengths for a region of interest (exons/introns)
#   #
#   # Args:
#   #   features:     The deseq2 object of interest (exons or introns)
#   #   region.type:  A string "exons" or "introns" for the name of file to get/make
#   #   dataset.name: Default is "gtex". The name of the dataset for which i want the lengths
#   #
#   # Returns:
#   #   A named vector with the lengths of the genes for each region of interest
#   #
#   
#   # Prepare the path to the folder with the gene length RDS
#   gene.lengths.folder <- paste(here("data-input"),
#                                "gene-lengths",
#                                sep = "/")
#   
#   # If the folder does not exist, create it
#   if(dir.exists(gene.lengths.folder) == FALSE) {
#     dir.create(gene.lengths.folder)
#   }
#   
#   # Create the path to the file
#   gene.lengths.path <- paste(gene.lengths.folder, 
#                              paste(dataset.name,
#                                    region.type,
#                                    "gene-lengths.rds",
#                                    sep = "-"),
#                              sep = "/")
#   # If the gene lengths file exist
#   if(file.exists(gene.lengths.path) == TRUE) {
#     
#     message("Loading gene length file for:", region.type)
#     
#     # Just load it    
#     gene.lengths <- readRDS(gene.lengths.path)
#     
#   } else {
#     
#     message("Generating gene length file for:", region.type, "")
#     
#     # Otherwise get the lengths
#     gene.lengths <- mclapply(rowRanges(features), function(gene) sum(gene@ranges@width))
#     
#     # And save them
#     saveRDS(gene.lengths, gene.lengths.path)
#   }
#   
#   return(gene.lengths)
# }

normalize.features.by.TPM <- function(features, region.type, ...) {
  
  #
  # Normalizes the features of interest (exons or introns) using the TPM method
  #
  # Args:
  #   features:     The deseq2 object of interest (exons or introns)
  #   region.type:  A string either "exons" or "introns" to get the corresponding gene lengths file
  #   ...:          Extra arguments for the name of gene length file
  #
  # Returns:
  #   A data.table with the normalized features of interest
  #
  
  # Formula taken from here
  # https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
  
  # Get the total length of exons or introns for each gene
  gene.lengths <- get.gene.info(features, region.type, ...)
  
  gene.lengths <- gene.lengths[names(rowRanges(features)) %in% names(gene.lengths)]
  
  # Get the counts for exons or introns and wrap them in a data.table
  features.datatable <- as.data.table(counts(features)) 
  
  # Turm the bases to kilobases
  gene.lengths.kb <- unlist(gene.lengths)/10^3 
  
  # Get the sample columns
  sample.columns <- grep("gene.id", colnames(features.datatable), value = TRUE, invert = TRUE)
  
  # Get the reads per kilobases
  features.datatable[, (sample.columns) := mclapply(.SD, function(x) x/gene.lengths.kb), .SDcols = sample.columns]
  
  # Get the transcript per million
  features.datatable[, (sample.columns) := mclapply(.SD, function(x) x/(sum(x)/10^6)), .SDcols = sample.columns]
  
  # Put the gene names as a column
  features.datatable[, gene.id := rownames(counts(features))]
  
  # Reorder the columns
  setcolorder(features.datatable, c("gene.id", sample.columns))
  
  return (features.datatable)
}

export.EI.ratios.with.annotation <- function(data, data.annotation, sample.columns,
                                             analysis.name, file.name = "EI-ratios-annotated.csv") {
  
  #
  # Takes the EI table, transposes it, merges it with the tissue annotation information 
  # and saves it.
  #
  # Args:
  #   data:           The data.table with the E/I ratios 
  #                   (columns are the sample names/rows are the gene names)
  #   data.annotation:The data.table with the annotation of the samples
  #   sample.columns: A vector with the column names of the samples
  #   analysis.name:  The string of the analysis.name (we use it as folder name)
  #   file.name:      The name of the exported file
  #
  # Returns:
  #   The saved EI-ratios-annotated.csv
  #
  
  # Traspose the data.table by gene name
  data.transposed <- dcast(melt(data[, .SD, .SDcols = c(sample.columns, "gene.id")],
                                id.vars = "gene.id"),
                           variable ~ gene.id)
  
  # Update the variable column name to sample.id
  setnames(data.transposed, "variable", "sample.id")
  
  # Merge the tissue annotation data.table with the E/I ratios data.table
  EI.ratios.with.annotation <- merge(data.transposed, data.annotation, by = "sample.id")
  
  # Constract the path for the csv file
  file.name.path <- here("data-output",
                         analysis.name,
                         file.name)
  
  message("Exporting E/I ratios with annotation csv in:", file.name.path ,"...")
  
  # And finally save the table
  fwrite(EI.ratios.with.annotation,
         file.name.path,
         sep = ",",
         quote = FALSE)
  
  message("Export complete!")
}
