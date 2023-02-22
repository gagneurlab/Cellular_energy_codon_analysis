# All functions exist in this file. If they are too many, they can be separated in
# several functions_XXX.R files.

# All functions exist in this file. If they are too many, they can be separated in
# several functions_XXX.R files.

# Return the memory to the OS
gc(verbose = FALSE,
   reset = TRUE)

# Functions scripts to load
functions.subfiles <- c("functions_load_data.R",
                        "functions_build.R"#,
                        #"functions_plots.R",
                        #"functions_model.R",
                        #"functions_sequences.R",
                        #"functions_analyse.R"
                        )

# Load the functions
invisible(lapply(functions.subfiles, source))
