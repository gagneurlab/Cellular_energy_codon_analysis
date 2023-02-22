check.parameters.validity <- function(analysis.parameters) {
  #
  # Reads the parameters and checks if they are valid, if not stops the analysis
  # and prints an appropriate message
  #
  # Requires:
  #
  # Load command:
  # 
  #
  # Args:
  #   analysis.parameters: A list the the parameters of the analysis
  #
  # Returns:
  #   Stops the analysis if an argument is not valid
  #
  
  # The supported organisms
  supported.organisms <- c("Homo sapiens",
                           "Mus musculus",
                           "Arabidopsis thaliana")
  
  # The supported genome versions
  supported.genome.versions <- c("hg19", "hg38",
                                 "mm9", "mm10",
                                 "TAIR10")
  
  # Check validity of arguments
  if(as.logical(analysis.parameters[["get.protein.coding"]] %in% c(TRUE, FALSE)) == FALSE) {
    
    stop("Invalid value for --get.protein.coding argument.\n",
         "Only accepted values are TRUE/FALSE.", call. = FALSE)
  }
  
  # Shoud we get all the samples from a tissue or just a portion?
  if(analysis.parameters[["samples.per.condition"]] != "all") {
    # Make sure that the samples.per.condition is a number and then typecast it
    # otherwise give stop the analysis
    if(suppressWarnings(is.na(as.integer(analysis.parameters[["samples.per.condition"]])) != TRUE &
                        as.numeric(analysis.parameters[["samples.per.condition"]]) %% 1  == 0 &
                        as.numeric(analysis.parameters[["samples.per.condition"]]) > 0)) {
      
      analysis.parameters[["samples.per.condition"]] <- as.numeric(analysis.parameters[["samples.per.condition"]])
      
      
    } else {
      stop("Invalid value for --samples.per.condition argument.\n",
           "Only accepted values are all or positive numbers greater than 0.", call. = FALSE)
    }  
  }
  
  # Check that the TPM argument is a number
  if(suppressWarnings(is.na(as.integer(analysis.parameters[["tpm.threshold"]])) != TRUE &
                      as.numeric(analysis.parameters[["tpm.threshold"]]) %% 1  == 0 &
                      as.numeric(analysis.parameters[["tpm.threshold"]]) >= 0)) {
    
    analysis.parameters[["tpm.threshold"]] <- as.numeric(analysis.parameters[["tpm.threshold"]])
    
    
  } else{
    stop("Invalid value for --tpm.threshold argument.\n",
         "Only accepted values are all or positive numbers greater or equal than/or 0.", call. = FALSE)
  }
  
  if(analysis.parameters[["tpm.table.path"]] != ""){
    # Subset the list of parameters for the file paths
    file.paths.subset <- analysis.parameters[names(analysis.parameters) %in% c("annotation.table.path",
                                                                               "tpm.table.path")]
  } else {
    # Subset the list of parameters for the file paths
    file.paths.subset <- analysis.parameters[names(analysis.parameters) %in% c("annotation.table.path")]
  }
  
  
  # Check that the path provided is a valid formated file path to a csv/txt/gz
  file.paths.are.valid <- grepl("^\\/?([[:alnum:]\\_\\-]+\\/)+([[:alnum:]\\_\\-\\.])+\\.[a-z]{2,3}$",
                                file.paths.subset, perl = TRUE)
  
  # Check that the path provided is a valid formated file path to a folder
  bam.files.directory.is.valid <- grepl("^\\/?([[:alnum:]_\\-]+\\/)*[[:alnum:]_\\-]+$",
                                        analysis.parameters[["bam.files.directory"]])
  
  if(all(file.paths.are.valid) == FALSE) {
    stop("Invalid value for ", paste("--", names(file.paths.subset)[!file.paths.are.valid],
                                     sep = "", collapse = ", ")," .\n",
         "Only accepted values are valid file paths", call. = FALSE)
  }
  
  # Is the bam file directory valid?
  if(bam.files.directory.is.valid != TRUE) {
    stop("Invalid value for --bam.files.directory .\n",
         "Only accepted values are valid folder paths", call. = FALSE)
  }
  
  # Is the organism provided in the supported organisms?
  if(any(analysis.parameters[["genome.version"]] %in% supported.genome.versions) == FALSE) {
    stop("Invalid value for --genome.version .\n",
         "Only supported genome versions are ", paste(supported.genome.versions, collapse = ", "),
         call. = FALSE)
    
  }
  
  # Is the genome version provided in the supported versions?
  if(any(analysis.parameters[["organism"]] %in% supported.organisms) == FALSE) {
    stop("Invalid value for --organism .\n",
         "Only supported organisms are ", paste(supported.organisms, collapse = ", "),
         call. = FALSE)
    
  }
}


check.command.line.arguments <- function(arguments = NULL,
                                         tolerate.no.arguments = TRUE,
                                         report.arguments.only = FALSE) {
  
  #
  # Checks if the correct arguments are supplied, in the correct order
  # and amount, and prints informative messages in case they are not
  #
  # Arguments:
  #   arguments:              Default is NULL. A vector with the argument names and 
  #                           the argument values
  #   tolerate.no.arguments:  Default is TRUE. If tolerate.no.arguments is set TRUE
  #                           I will not stop the analysis when no arguments are provided
  #   report.arguments.only:  Default is FALSE. Helps getting the default arguments
  #
  # Return:
  #   Returns a list with the correct arguments for the analysis or the default
  #   (if some of them are not provided) or throws an error message and stops the analysis
  #
  
  # Read the default argument table
  default.argument.table <- data.table::fread(here("data-input/default-analysis-parameters.csv"))
  
  # Get the available arguments
  available.arguments <- default.argument.table[, parameter]
  
  # And their default values
  default.command.line.argument.values <- default.argument.table[, value]
  
  # # Define the available arguments
  # available.arguments <- c("dataset.type",
  #                             "analysis.name",
  #                             "get.protein.coding",
  #                             "is.paired.end",
  #                             "samples.file.name",
  #                             "samples.per.condition",
  #                             "condition.column.major",
  #                             "condition.column.minor",
  #                             "annotation.table.path",
  #                             "conditions.major",
  #                             "conditions.minor",
  #                             "sample.column",
  #                             "bam.files.directory",
  #                             "samples.to.discard.pattern",
  #                             "tpm.threshold",
  #                             "tpm.table.path",
  #                             "threads",
  #                             "organism",
  #                             "genome.version",
  #                             "gencode.version")
  # 
  # default.command.line.argument.values <- c("gtex",
  #                                           "mrna-half-life-analysis",
  #                                           TRUE,
  #                                           TRUE,
  #                                           "samples.file.name",
  #                                           "all",
  #                                           "conditions",
  #                                           "",
  #                                           "/annotation/table/path",
  #                                           "conditions.column",
  #                                           "",
  #                                           "RNA_ID",
  #                                           "/bam/files/directory",
  #                                           "",
  #                                           "1",
  #                                           "",
  #                                           "8",
  #                                           "Homo sapiens",
  #                                           "hg19",
  #                                           "19")
  
  # In case I want only to report the arguments
  if(report.arguments.only == TRUE) {
    return(available.arguments) 
  }
  
  # Find the argument names
  argument.names <- grep("--", arguments, value = TRUE)
  
  # And the argument name positions
  argument.names.positions <- grep("--", arguments)
  
  # Get the argument values
  argument.values <- grep("--", arguments, value = TRUE, invert = TRUE)
  
  # And the argument values positions
  argument.values.positions <- grep("--", arguments, invert = TRUE)
  
  # Now if the argument names are in the correct positions
  # They should occupy the even positions (1, 3, etc...)
  argument.names.correct.positions <- argument.names.positions[argument.names.positions %% 2 == 1]
  
  # Are they in the correct positions?
  argument.names.are.in.correct.positions <-  all(argument.names.positions %% 2 == 1)
  
  # Also if the argument values are in the correct positions
  # They should occupy the odd positions (2, 4, etc...)
  argument.values.correct.positions <- argument.values.positions[argument.values.positions %% 2 == 0]
  
  # Are they in the correct positions?
  argument.values.are.in.correct.positions <-  all(argument.values.positions %% 2 == 0)
  
  # First we will check that the arguments used are the ones that I should expect
  arguments.are.valid <- all(argument.names %in% paste0("--", available.arguments))
  
  # Is the argument vector empty or not?
  if(is.null(arguments) == TRUE | length(arguments) == 0) {
    
    # Should I tolerate if the user gave no arguments?
    # If yes perharps he will provide an analysis file
    if(tolerate.no.arguments == TRUE) {
      message("No analysis parameters provided by command line.")
      message("Will attemp to use the analysis-parameters.csv...")
      
      return(TRUE)
    } else {
      
      # In case we do not tolerate give an error message
      error.message.1 <- paste0("No command line arguments or analysis-parameters.csv input file provided.\n",
                                "Please provide arguments for the analysis.")
      
      error.message.2 <- c("Available parameters are the following:\n",
                           paste0("\t", "--", available.arguments,"\n"))
      
      # Finallize the error message
      error.message <- c(error.message.1,
                         error.message.2)
      
      stop(error.message, call. = FALSE)  
    }
    
  } else if(is.null(arguments) == FALSE) {
      
    if(length(arguments) != 0) {
      
          message("Command line arguments supplied. Checking validity...")
      
      } else if(tolerate.no.arguments == TRUE & length(arguments) == 0){
        
        message("No analysis parameters provided by command line.")
        message("Will attemp to use the analysis-parameters.csv...")
      } else {
        
        stop("Uncaught exception in check.command.line.arguments")
      
      }
    
  }
  
  # If the arguments are not valid
  if(arguments.are.valid == FALSE) {
    
    # Find the invalid ones
    invalid.argument <- argument.names[!argument.names %in% paste0("--", available.arguments)]
    
    # And inform the user with the appropriate messages
    error.message.1 <- paste0("No such valid argument(s) as: ",
                              paste(invalid.argument, collapse = ", "),
                              ".\n")
    
    error.message.2 <- c("Available arguments are the following:\n",
                         paste("\t", "--", available.arguments,"\n"),
                         ifelse(length(invalid.argument) == 1,
                                "Is it an argument value?\n",
                                "Are they argument values?\n"))
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
  }
  
  # Now get the number of duplicates
  number.of.duplicates <- sum(duplicated(argument.names))
  
  # If I have duplicates, print the approriate message
  if(number.of.duplicates == 1) {
    
    # Prepare the error messages
    error.message.1 <- paste0("Duplicated argument name: ",
                              paste(argument.names[duplicated(argument.names)], collapse = ", "),
                              " is duplicated.\n")
    
    error.message.2 <- c("Available arguments are the following:\n",
                         paste0("\t", "--", available.arguments,"\n"))
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
    
  } else if(number.of.duplicates > 1) {
    
    # Prepare the error messages
    error.message.1 <- paste0("Duplicated argument names: ",
                              paste(argument.names[duplicated(argument.names)], collapse = ", "),
                              " are duplicated.\n")
    
    error.message.2 <- c("Available arguments are the following:\n",
                         paste0("\t", "--", available.arguments,"\n"))
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
  }
  
  # Finally lets see if the arguments are in wrong positions
  # If the number of arguments names is equal to the argument values
  # but they are positioned in the wrong positions
  # Print the approriate message
  # Otherwise everything should be fine
  if(length(argument.names) == length(argument.values) & 
     (argument.names.are.in.correct.positions == FALSE |
      argument.values.are.in.correct.positions == FALSE)) {
    
    error.message.1 <- "Arguments are in incorrect order.\n"
    
    error.message.2 <- c("Order should be ",
                         paste0("--", available.arguments[1:2], " arg.value", collapse = ", "),
                         " etc.\n")
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
    
  }else if(length(argument.names) == length(argument.values) & 
           argument.names.are.in.correct.positions == TRUE &
           argument.values.are.in.correct.positions == TRUE) {
    
    # If the argument names are equal to the argument values
    # and they are positioned one by one (one name one value)
    # Everything is fine
    # message("All good")
    
    # Prepare the user parameters list
    user.parameters <- as.list(argument.values)
    
    # Add the names
    names(user.parameters) <- gsub("--", "", argument.names)
    
    # Prepare default parameters list
    default.parameters <- as.list(default.command.line.argument.values)
    
    names(default.parameters) <- available.arguments
    
    # Get not user provided arguments
    not.user.provided.arguments <- default.parameters[!names(default.parameters) %in%
                                                        names(user.parameters)]
    
    # Concatenate the user arguments with the default arguments
    merged.arguments <- c(user.parameters, not.user.provided.arguments)
    
    # And put the in the correct order
    merged.arguments.reordered <- merged.arguments[order(match(names(merged.arguments),
                                                             available.arguments))]
    
    if(length(user.parameters) == 1 &
       names(user.parameters)[1]  == "analysis.name") {
      message("--analysis.name provided. Will attempt to read the provided folder...")
      return(user.parameters)
    } else {
      return(merged.arguments.reordered)
    }
    
  } 
  
  # If nothing is caught till now
  # Either ther argument names are more than the argument values
  # or the other way around
  if(length(argument.names) < length(argument.values)) {
    
    # If the argument values are more than the argument names
    error.message.1 <- "Incorrect number of arguments. More argument values than argument names.\n"
    
    # Find the wrong arguments
    wrong.arguments <- arguments[argument.values.positions[argument.values.positions %% 2 == 1]]
    
    # If the arguments are wrong maybe is a typo
    typos <- available.arguments[which(available.arguments %in% arguments)]
    
    # Prepare the error message 2
    error.message.2 <- ""
    
    # Inform the user for the typos
    if(length(typos) > 0) {
      
      # For printing purposes get only a slice from the typos
      if(length(typos) > 1){
        
        typos.slice <- typos[1:2]
        
      } else {
        
        typos.slice <- typos
        
      }
      
      error.message.2 <- c(paste0(typos, collapse = ", "),
                           " seem to belong in the valid argument names though.\n",
                           "Did you forget to add '--' before these argument(s) e.g. ",
                           paste0("--", typos.slice, " arg.value", collapse = ", "),
                           "?\n")
    }
    
    # And the rest of the wrong arguments should be invalid
    invalid.arguments <- setdiff(wrong.arguments, typos)
    
    if(length(invalid.arguments) > 0) {
      
      error.message.2 <- c(error.message.2,
                           paste0("No such argument(s) as: ",
                                  paste(invalid.arguments, collapse = ", "),
                                  "."),
                           "\n")
    }
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
    
    return(FALSE)
    
  } else if(length(argument.names) > length(argument.values)) {
    
    error.message.1 <- "Incorrect number of arguments. More argument names than argument values.\n"
    
    # The argument values should be shifted by 1
    solo.arguments.positions <- setdiff(argument.names.positions, argument.values.positions-1)
    
    # Get the wrong arguments
    wrong.arguments <- arguments[solo.arguments.positions]
    
    # For printing purposes get only a slice from the wrong arguments
    if(length(wrong.arguments) > 1) {
      
      wrong.arguments.slice <- wrong.arguments[1:2]
      
    } else {
      
      wrong.arguments.slice <- wrong.arguments
    }
    
    error.message.2 <- c(paste0(wrong.arguments, collapse = ", "),
                         ifelse(length(wrong.arguments) == 1,
                                " does not",
                                " do not"),
                         " have argument values.\n",
                         "Argument(s) should be e.g. ",
                         paste0(wrong.arguments.slice, " arg.value", collapse = ", "),
                         " etc.\n")
    
    # Finallize the error message
    error.message <- c(error.message.1,
                       error.message.2)
    
    stop(error.message, call. = FALSE)
  } 
}

check.and.install.packages <- function(packages, origin = "CRAN") {
  #
  # Checks to see if desired packages are installed and if they are not, it installs them.
  #
  # Args:
  #   package:  A vector of the desired packages
  #   origin:   Default is CRAN. Packages can be from CRAN, GitHub or Bioconductor
  #
  # Returns:
  #   TRUE by default
  #
  
  switch (origin,
          
          "CRAN" = {
            new.packages <- packages[ !(packages %in% installed.packages()[, "Package"])]
          },
          "GitHub" = {
            # For github packages each package is in the form repositoryName/packageName,
            # so in order to see if the package exists, I split the input on the slash
            # and get every second element
            splited.packages <- unlist(strsplit(packages, "/"))[c(FALSE, TRUE)]
            
            new.packages <- packages[ !(splited.packages %in% installed.packages()[, "Package"])]
          },
          "Bioconductor" = {
            new.packages <- packages[ !(packages %in% installed.packages()[, "Package"])]
          })
  
  if (length(new.packages) > 0) {
    switch (origin,
            "CRAN" = {
              install.packages(new.packages, dependencies = TRUE)
            },
            "GitHub" = { 
              install_github(new.packages)
            },
            "Bioconductor" = {
              install(new.packages)
            }) 
  }  
  return (TRUE)
}
