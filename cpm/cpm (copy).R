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
output_path <- getOptionValue("--output_path")
input_file1 <- getOptionValue("--galaxy_input1")
#input_file2 <- getOptionValue("--galaxy_input2")
output_file1 <- getOptionValue("--galaxy_output1")
#output_file2 <- getOptionValue("--galaxy_output2")
param1 <- getOptionValue("--galaxy_param1")
param2 <- getOptionValue("--galaxy_param2")
param3 <- getOptionValue("--galaxy_param3")
param4 <- getOptionValue("--galaxy_param4")
param5 <- getOptionValue("--galaxy_param5")

# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)


if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}


library(MASS)
library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(tidyr)


# Processing input files
file_list <- strsplit(input_file1, ",")[[1]]

# Remove specific files if param3 is not empty
if (param3 != "") {
  rmrun_list <- as.list(strsplit(param3, ",")[[1]])
  for (name in rmrun_list) {
    rmrun_file <- grep(pattern = paste0(name, "_count.txt$"), file_list)
    if (length(rmrun_file) > 0) {
      file_list <- file_list[-rmrun_file]
    }
  }
}

# Read gene names from the first file
gene_name <- read.delim(file_list[[1]], header = FALSE)
dataframe_total_count <- data.frame(Geneid = gene_name[, 1])

# Combine counts from all files
for (file in file_list) {
  count <- read.delim(file, header = FALSE, comment.char = "#", quote = "")
  count[1, ] <- sub(".*\\/.*\\/", "", count[1, ])
  dataframe_total_count <- cbind(dataframe_total_count, count[, 2])
}


# Write the combined count table
write.table(dataframe_total_count, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Paths for input and coldata
count_filtered_path <- output_file1
coldata_path <- param5

# Load data
count_filtered_df <- read_delim(count_filtered_path, col_names = FALSE)

# Load metadata
coldata_df <- read_delim(coldata_path)

# Update column names
colnames(count_filtered_df)[2:ncol(count_filtered_df)] <- coldata_df$project

# Detect unique conditions and create combinations
conditions <- unique(coldata_df$condition)
condition_combinations <- combn(conditions, 2, simplify = FALSE)

#counter <- 0

#file_list_output <- list()
#output_file2 <- file_list_output
#output_directory <- "/home/audrey/Documents/Bulk_RNAseq_analyses/Result_galaxy"


# Process each condition combination
for (combo in condition_combinations) {
  #counter <- counter + 1

  condition1 <- combo[1]
  condition2 <- combo[2]
  
  samples_condition1 <- coldata_df %>% filter(condition == condition1) %>% pull(project)
  samples_condition2 <- coldata_df %>% filter(condition == condition2) %>% pull(project)
  
  comparison_df <- count_filtered_df %>%
    select(X1, all_of(samples_condition1), all_of(samples_condition2)) %>%
    rename(Geneid = X1)

  #file_name <- paste0(condition1, ".tabular")
  name <- paste0(condition1, "_vs_", condition2)
  file_path <- file.path(output_path, paste0(name, ".tabular")) 
  write.table(comparison_df, file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  #file_list_output[[counter]] <- file_path
}




