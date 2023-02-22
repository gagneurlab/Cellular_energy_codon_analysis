# Data wrangling with dplyr/tidyr/etc. All the magic happens here. Newly created
# files in workspace should be displayed. A build boolean variable can be used for
# data reload from RDS for faster data reload.

# Clear enviroment and only keep functions and global/project variables
rm(list = grep(paste(c("^global.variables",
                       "^project.variables",
                       lsf.str()),
                     collapse = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Return the memory to the OS
gc(verbose = FALSE,
   reset = TRUE)

# Get the sample table
samples.table <- global.variables[["samples.table"]] 

# Are the samples pair ended?
is.paired.end <- global.variables[["is.paired.end"]]

# Get the TPM threshold
tpm.threshold <- global.variables[["tpm.threshold"]]

# Set it to number
global.variables[["tpm.threshold"]] <-  as.numeric(global.variables[["tpm.threshold"]])

# How many threads to use for multithread processes
threads.to.use <- global.variables[["threads"]]

# Get the organism
organism <- global.variables[["organism"]]

# Get the exons, introns and raw exons granges
exons <- global.variables[["exons.granges"]]
introns <- global.variables[["introns.granges"]]
raw.exons <- global.variables[["raw.exons.granges"]]

# Get the analysis name
analysis.name <- global.variables[["analysis.name"]]

# Get the annotation table path
annotation.table.path <- global.variables[["annotation.table.path"]]

# Get the conditions to compare
conditions.to.compare <- global.variables[["conditions"]]

# Get the path of the TPM table
tpm.table.path <- global.variables[["tpm.table.path"]]

# Get the column name of the column containing the major conditions
condition.column.major <- global.variables[["condition.column.major"]]

# Get the column name of the column containing the minor conditions
condition.column.minor <- global.variables[["condition.column.minor"]]

# Find the sample ids
srr.id.paths <- samples.table[, RNA_BAM_FILE]

# Get the genome version
genome.version <- global.variables[["genome.version"]]

# Get the gene annotation version
genes.annotation.version <- global.variables[["gencode.version"]]

# Prepare the bam files paths
samples.to.condition <- samples.table[, .SD, .SDcols = c("RNA_ID",
                                                         condition.column.major,
                                                          "RNA_BAM_FILE")]

# Add the SRR sample ids as names on the path vector
names(srr.id.paths) <- samples.table[, RNA_ID]

# Do we have the exons/introns objects?
# Prepare the RDS file name
rds.file.name <- paste0("exons",
                        "-counts-",
                        analysis.name,
                        ".rds")

# Prepare the RDS file path
rds.file.path <- here("data-output",
                      analysis.name,
                      "counts-rds",
                      rds.file.name)

message("Checking if samples are paired or single ended...")
# 
# # Check if the bams are paired end
# bams.are.paired.ended <- mcmapply(testPairedEndBam,
#                                   srr.id.paths,
#                                   mc.cores = threads.to.use)
# 
# # Check if everything is the same
# if(sum(bams.are.paired.ended) == 0 |
#    sum(bams.are.paired.ended) == length(bams.are.paired.ended)) {
# 
#   # If no sample is Paired-end set the flag to FALSE
#   if(sum(bams.are.paired.ended) == 0) {
#     is.paired.end <-  FALSE
# 
#     message("BAM files are single-end")
#   } else {
#     is.paired.end <-  TRUE
#     message("BAM files are pair-end")
# 
#   }
# } else {
#   stop("Mixed paired end and single end samples!")
# }

# If RDS object exists skip the BamFileList creation
if(file.exists(rds.file.path) == FALSE) {
  
  message(Sys.time(), " Creating the BamFilesList")
  
  # Make a bam files list
  bam.file.list <- BamFileList(srr.id.paths,
                               asMates = is.paired.end,
                               yieldSize = 2e6)
  
  # Set the bam files seq style to Ensembl style (1, 2, ..., X, Mt, Pt) 
  # seqlevelsStyle(bam.file.list) <- "Ensembl"
  
  message(Sys.time(), " Setting sequence style to Ensembl")
  
  # Set introns seq style to Ensembl style (1, 2, ..., X, Mt, Pt)
  seqlevelsStyle(exons) <- seqlevelsStyle(bam.file.list[1])
  
  # Set introns seq style to Ensembl style (1, 2, ..., X, Mt, Pt) 
  seqlevelsStyle(introns) <- seqlevelsStyle(bam.file.list[1])
}

# Get the count table for exons
exons.count.table <- get.counts.table(features = exons,
                                      analysis.name = analysis.name,
                                      table.type = "exons",
                                      samples.to.condition = samples.to.condition,
                                      condition.column = condition.column.major,
                                      bam.file.list = bam.file.list,
                                      is.paired.end = is.paired.end,
                                      count.fragments = FALSE,
                                      # overwrite = TRUE,
                                      number.of.threads = threads.to.use)

stop("end test")
# Get the count table for introns
introns.count.table <- get.counts.table(features = introns,
                                        analysis.name = analysis.name,
                                        table.type = "introns",
                                        samples.to.condition = samples.to.condition,
                                        condition.column = condition.column.major,
                                        bam.file.list = bam.file.list,
                                        is.paired.end = is.paired.end,
                                        count.fragments = FALSE,
                                        # overwrite = TRUE,
                                        number.of.threads = threads.to.use)

# Count the raw exons
raw.exons.counts <- get.counts.table(features = raw.exons,
                                     analysis.name = analysis.name,
                                     table.type = "raw-exons",
                                     samples.to.condition = samples.to.condition,
                                     condition.column = condition.column.major,
                                     bam.file.list = bam.file.list,
                                     is.paired.end = is.paired.end,
                                     count.fragments = FALSE,
                                     # overwrite = TRUE,
                                     number.of.threads = threads.to.use)

# If there is no TPM file provided, get the raw exons counts as well, generate the TPM table,
# And renew the path
if(tpm.table.path == "") {

  # Create the TPM table per sample and get the path
  tpm.table.path <- generate.tpms(raw.exon.counts = raw.exon.counts,
                                  analysis.name = analysis.name,
                                  genome.version = genome.version,
                                  # overwrite = TRUE,
                                  genes.annotation.version = genes.annotation.version)

  # Update the path of the TPM table
  global.variables[["tpm.table.path"]] <- tpm.table.path
}

# Store raw counts of exons
raw.exon.counts <- store.raw.counts(exons.count.table, "exons", analysis.name,
                                    # overwrite = TRUE
                                    )

# Store raw counts of introns
raw.intron.counts <-store.raw.counts(introns.count.table, "introns", analysis.name,
                                     # overwrite = TRUE
                                      )
# Store raw counts of exons
invisible(store.raw.counts(features = raw.exons.counts,
                           feature.type = "raw-exons",
                           analysis.name = analysis.name,
                           # protein.coding.only = protein.coding.only,
                           # overwrite = TRUE,
                           return = FALSE))

# Fix the levels of the conditions factors by removing weird charachters
exons.count.table <- fix.SummarizedExperiment.condition.names(exons.count.table,
                                                               condition.column.major)

# Fix the levels of the conditions factors by removing weird charachters
introns.count.table <- fix.SummarizedExperiment.condition.names(introns.count.table,
                                                                 condition.column.major)


# Fix the levels of the conditions factors by removing weird charachters
raw.exons.counts <- fix.SummarizedExperiment.condition.names(raw.exons.counts,
                                                             condition.column.major)

# Do I have multiple tissues or only one?
if(length(conditions.to.compare) > 1) {

  # Wrap exons as a DESeq object
  exons.deseq2 <- DESeqDataSet(exons.count.table, as.formula(paste("", condition.column.major, sep = "~")))

  # Wrap introns as a DESeq object
  introns.deseq2 <- DESeqDataSet(introns.count.table, as.formula(paste("", condition.column.major, sep = "~")))
} else {

  # Wrap exons as a DESeq object
  exons.deseq2 <- DESeqDataSet(exons.count.table, ~ 1)

  # Wrap introns as a DESeq object
  introns.deseq2 <- DESeqDataSet(introns.count.table, ~ 1)
}

# Remove redundant objects
rm(exons.count.table, introns.count.table)

# Get the exons/introns ratios
EI.ratios.centered <- combine.exons.introns(exons.deseq2 = exons.deseq2,
                                            introns.deseq2 = introns.deseq2,
                                            analysis.name = analysis.name,
                                            annotation.table.path = annotation.table.path,
                                            condition.column = condition.column.major,
                                            tpm.table.path = tpm.table.path,
                                            number.of.threads = threads.to.use,
                                            # overwrite = TRUE,
                                            # tpm.threshold = 0.5,
                                            genome.version = genome.version)

# Get sample columns
sample.columns <- grep("gene.id|gene.average",
                       colnames(EI.ratios.centered), value = TRUE, invert = TRUE)

# Make a condition annotation data.table
condition.annotations <- mapply(get.samples.annotation,
                                condition.column = c(condition.column.major,
                                                     condition.column.minor),
                                MoreArgs = list(features.deseq2 = colData(exons.deseq2),
                                               analysis.name = analysis.name,
                                               # overwrite = TRUE,
                                               annotation.table.path = annotation.table.path))

# Prepared the mean E/I ratio per condition (major/minor)
EI.ratios.per.conditions <- mapply(export.average.EI.by.condition,
                                   condition.annotation = condition.annotations,
                                   condition.column = c(condition.column.major,
                                                        condition.column.minor),
                                   MoreArgs = list(EI.ratios.centered = EI.ratios.centered,
                                                   sample.columns = sample.columns,
                                                   # overwrite = TRUE,
                                                   analysis.name = analysis.name))

# Load the non centered EI ratios
EI.ratios.non.centered <- fread(here("data-output",
                                     analysis.name,
                                     "EI-ratios-masked.csv"))

# Prepared the mean E/I ratio per condition for not centered data (major/minor)
EI.ratios.per.conditions.non.centered <- mapply(export.average.EI.by.condition,
                                               condition.annotation = condition.annotations,
                                               condition.column = c(condition.column.major,
                                                                    condition.column.minor),
                                               MoreArgs = list(EI.ratios.centered = EI.ratios.non.centered,
                                                               sample.columns = sample.columns,
                                                               analysis.name = analysis.name,
                                                               # overwrite = TRUE,
                                                               file.name.suffix = "EI-ratios-non-centered-per"))


# Add the sample.columns as global.variable
global.variables[["sample.columns"]] <- sample.columns

# Add the EI ratios per sample table as global.variable
global.variables[["EI.ratios.centered"]] <- EI.ratios.centered

# Add the EI ratios per sample non centered table as global.variable
global.variables[["EI.ratios.non.centered"]] <- EI.ratios.non.centered

# Add the condition.annotation table as global.variable
global.variables[["condition.annotation.major"]] <- condition.annotations[[1]]

# Add the EI.per condition table as global.variable
global.variables[["EI.ratios.per.condition.major"]] <- EI.ratios.per.conditions[[1]]

# Add the EI.per condition non centered table as global.variable
global.variables[["EI.ratios.per.condition.non.centered.major"]] <- EI.ratios.per.conditions.non.centered[[1]]

# Add the minor condition annotation if minor condition column was provided
if(is.null(condition.column.minor) == FALSE) {
  
  # Add the condition.annotation table as global.variable
  global.variables[["condition.annotation.minor"]] <- condition.annotations[[2]]
  
  # Add the EI.per tissue table as global.variable
  global.variables[["EI.ratios.per.condition.minor"]] <- EI.ratios.per.conditions[[2]]
  
  # Add the EI.per condition non centered table as global.variable
  global.variables[["EI.ratios.per.condition.non.centered.minor"]] <- EI.ratios.per.conditions.non.centered[[2]]
}
  
cat("========== End of build.R ==========\n")

