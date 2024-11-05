# How to execute this tool
# $Rscript my_r_tool.R --input input.txt --output output.txt

# Send R errors to stderr
options(show.error.messages = FALSE, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, FALSE)
})

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
input_file1 <- getOptionValue("--galaxy_input1")
#input_file2 <- getOptionValue("--galaxy_input2")
output_file1 <- getOptionValue("--galaxy_output1")
output_file2 <- getOptionValue("--galaxy_output2")
output_file3 <- getOptionValue("--galaxy_output3")
param1 <- getOptionValue("--galaxy_param1")
param2 <- getOptionValue("--galaxy_param2")
param3 <- getOptionValue("--galaxy_param3")
param4 <- getOptionValue("--galaxy_param4")
param5 <- getOptionValue("--galaxy_param5")

# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)



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

# Write the combined count table with column names as the first row
noms_colonnes <- names(dataframe_total_count)
write.table(t(noms_colonnes), file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Append the data to the same file
write.table(dataframe_total_count, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

# Paths for input and coldata
count_filtered_path <- output_file1
coldata_path <- param5

# Load data
count_filtered_df <- read_delim(count_filtered_path, col_names = FALSE)
coldata_df <- read_delim(coldata_path)

# Update column names
colnames(count_filtered_df)[2:ncol(count_filtered_df)] <- coldata_df$project

# Detect unique conditions and create combinations
conditions <- unique(coldata_df$condition)
condition_combinations <- combn(conditions, 2, simplify = FALSE)


# Process each condition combination
for (combo in condition_combinations) {
  condition1 <- combo[1]
  condition2 <- combo[2]
  
  samples_condition1 <- coldata_df %>% filter(condition == condition1) %>% pull(project)
  samples_condition2 <- coldata_df %>% filter(condition == condition2) %>% pull(project)
  
  comparison_df <- count_filtered_df %>%
    select(X1, all_of(samples_condition1), all_of(samples_condition2)) %>%
    rename(Geneid = X1)
  
  file_name <- paste0(condition1, "_vs_", condition2, ".tabular")

  # Write column names as the first row
  comparison_col_names <- names(comparison_df)
  write.table(t(comparison_col_names), file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Append the data to the same file
  write_tsv(comparison_df, file_name, append = TRUE)

  # Matrix transformation for CPM calculation
  mtx_total_count <- as.matrix(comparison_df[-1, -1])
  class(mtx_total_count) <- "numeric"
  cpm_values <- cpm(mtx_total_count)
  df_cpm <- as.data.frame(cpm_values)
  df_cpm <- cbind(Geneid = comparison_df[-1, 1], df_cpm)

  df_cpm[] <- lapply(df_cpm, function(x) if(is.factor(x)) as.numeric(as.character(x)) else x)

  # Filter CPM values
  df_cpm_filter <- df_cpm[apply(df_cpm[, -1], 1, function(row) sum(row >= param1) >= param2), ]
  df_cpm_filter <- rbind(comparison_df[1, ], df_cpm_filter)
  write.table(df_cpm_filter, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Filter counts based on filtered CPM
  filtered_genes <- df_cpm_filter$Geneid[-1]
  dataframe_filtered_count <- comparison_df[comparison_df$Geneid %in% filtered_genes, ]
  write.table(dataframe_filtered_count, file = output_file3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
