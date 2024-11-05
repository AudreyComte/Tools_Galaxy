# #################################################################
# 
#                               DESeq2
#                          
# #################################################################


# How to execute this tool
# $Rscript my_r_tool.R --input input.txt --output output.txt

# Send R errors to stderr
#options(show.error.messages = FALSE, error = function() {
  #cat(geterrmessage(), file = stderr())
  #q("no", 1, FALSE)
#})

# Avoid crashing Galaxy with a UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to get the value of an option from command line arguments
getOptionValue <- function(option) {
  idx <- match(option, args)
  if (!is.na(idx) && (idx + 1) <= length(args)) {
    return(args[idx + 1])
  }
  return(NULL)
}


# Get options
output_path2 <- getOptionValue("--output_path2")
output_path1 <- getOptionValue("--output_path1")


input <- getOptionValue("--galaxy_input")

output_names_file <- getOptionValue("--galaxy_output_names_files")
file.create(output_names_file)

output_names_ref <- getOptionValue("--galaxy_output_names_ref")
file.create(output_names_ref)

output_names_mutant <- getOptionValue("--galaxy_output_names_mutant")
file.create(output_names_mutant)

thread_param <- getOptionValue("--galaxy_param_thread")
param2_coldata <- getOptionValue("--galaxy_param2")
#param3_rmproj_list <- getOptionValue("--galaxy_param3")

# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)

if (!dir.exists(output_path2)) {
  dir.create(output_path2, recursive = TRUE)
}

if (!dir.exists(output_path1)) {
  dir.create(output_path1, recursive = TRUE)
}



########################################################################################################################################################################################################################################################################################################################################################################################################


library(DESeq2)
library(readr)

#parallel <- FALSE
#if (thread_param > 1) {
  #library("BiocParallel")
  #setup parallelization
  #register(MulticoreParam(thread_param))
  #parallel <- TRUE
#}


# Split the input string into a list of file paths
file_list <- strsplit(input, ",")[[1]]

# Read colData
coldata_read <- read_tsv(param2_coldata)

# Function to process a single file
process_file <- function(file_name) {
  # Read and convert the counts matrix
  cts <- as.matrix(read.table(file_name, header = TRUE, row.names = 1))
  
  # Filter rows in colData where the first column matches the headers in the counts matrix
  matching_rows <- coldata_read[coldata_read[[1]] %in% colnames(cts), ]
  
  # Create colData DataFrame from filtered rows
  coldata_filtered <- data.frame(matching_rows)
  
  # Set row names in colData
  rownames(coldata_filtered) <- coldata_filtered$project
  
  # Remove the 'project' column if it is no longer needed
  coldata <- coldata_filtered[, -which(names(coldata_filtered) == "project")]
  
  # Convert relevant columns to factors
  coldata$condition <- factor(coldata$condition)
  coldata$type <- factor(coldata$type)
  
  # Create the DESeqDataSet object
  if (ncol(cts) == nrow(coldata)) {
    dds <- DESeqDataSetFromMatrix(countData = cts, 
                                  colData = coldata, 
                                  design = ~ condition)
  } else {
    stop("Sample names don't match in both files: ", file_name, ".")
  }
  
  # Retrieve all levels of the condition
  all_levels <- levels(dds$condition)
  
  # Function to process each reference level
  process_level <- function(ref_level) {
    
    # Relevel the condition factor
    dds$condition <- relevel(dds$condition, ref = ref_level)
    
    # Get the other condition name
    mutant_level <- setdiff(all_levels, ref_level)[1]  # Assuming there is only one other level
    
    # DESeq: Normalization and preprocessing
    dds_result <- DESeq(dds)

    # Save the DESeq object
    name <- paste0(ref_level, "_vs_", mutant_level)
    output_file2 <- file.path(output_path2, paste0(name, ".rds"))
  
    #output_file2 <- file.path(output_path2, paste0("all_rds_", ref_level, "_vs_", mutant_level, ".rds"))
    saveRDS(dds_result, file = output_file2)
    
    # Append the name of the output file to the text file
    write(name, file = output_names_file, append = TRUE)
    write(ref_level, file = output_names_ref, append = TRUE)
    write(mutant_level, file = output_names_mutant, append = TRUE)

    # Save the normalized data matrix
    normalized_counts <- counts(dds_result, normalized = TRUE)

    name <- paste0(ref_level, "_vs_", mutant_level)
    output_file1 <- file.path(output_path1, paste0(name, ".tsv"))

    #output_file1 <- file.path(output_path1, paste0("normalized_counts_", ref_level, "_vs_", mutant_level, ".tsv"))
    write.table(normalized_counts, file = output_file1, sep = "\t", quote = FALSE, col.names = NA)
  }
  
  # Apply the function to each level without an explicit loop
  lapply(all_levels, process_level)
}

# Apply the function to each file in file_list
lapply(file_list, process_file)
