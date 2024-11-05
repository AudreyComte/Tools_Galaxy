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

# Directory output
output_path <- getOptionValue("--output_path")
output_cpm_path <- getOptionValue("--output_cpm_path")
output_filter_count_path <- getOptionValue("--output_filter_count_path")

# Input
input_file1 <- getOptionValue("--galaxy_input1")

# Output
output_file1 <- getOptionValue("--galaxy_output1")

# Parameter
param1 <- getOptionValue("--galaxy_param1")
param2 <- getOptionValue("--galaxy_param2")
param3 <- getOptionValue("--galaxy_param3")
param4 <- getOptionValue("--galaxy_param4")
param5 <- getOptionValue("--galaxy_param5")


# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)


# Creation of the output_path directory
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Creation of the output_path directory
if (!dir.exists(output_cpm_path)) {
  dir.create(output_cpm_path, recursive = TRUE)
}

# Creation of the output_path directory
if (!dir.exists(output_filter_count_path)) {
  dir.create(output_filter_count_path, recursive = TRUE)
}


library(MASS)
library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(tidyr)

# Processing input files
file_list <- strsplit(input_file1, ",")[[1]]

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

# Load data
output_count <- output_file1
count_filtered_df <- read_delim(output_count, col_names = FALSE)

# Load metadata
coldata <- param5
coldata_df <- read_delim(coldata)

# Copy the first row to use as column name
column_names <- as.character(unlist(count_filtered_df[1, 2:ncol(count_filtered_df)]))
colnames(count_filtered_df)[2:ncol(count_filtered_df)] <- column_names

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

  name <- paste0(condition1, "_vs_", condition2)
  output_count_comparaison <- file.path(output_path, paste0(name, ".tabular")) 
  write.table(comparison_df, output_count_comparaison, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Matrix transformation for cpm calculation
  mtx_total_count <- as.matrix(comparison_df[-1,-1])
  class(mtx_total_count) <- "numeric"
  cpm <- cpm(mtx_total_count)
  df_cpm <- as.data.frame(cpm)
  df_cpm <- cbind(Gene_id=comparison_df[-1,1], df_cpm)

  df_cpm[] <- lapply(df_cpm, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  # Dataframe with the cpm filtrated : It keep the genes with cpm value >= x in at least n samples
  # if n = 3 in a replicate keep the gene even if one experiment have 0 in each replicate
  df_cpm_filter <- data.frame()
  for (j in 1:nrow(df_cpm)) {
    flag = 0
    for (k in 2:ncol(df_cpm)) {
      if(df_cpm[j,k] >= param1){
        flag <- flag+1
      }
    }
    
    if(flag >= param2){
        df_cpm_filter <- rbind(df_cpm_filter,df_cpm[j,])
    }
  }

  df_cpm_filter <- rbind(comparison_df[1,],df_cpm_filter)
  name_cpm_filter <- paste0("output_cpm_", condition1, "_vs_", condition2)
  output_cpm_filter_comparaison <- file.path(output_cpm_path, paste0(name_cpm_filter, ".tabular")) 
  write.table(df_cpm_filter, output_cpm_filter_comparaison, sep=" ", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Dataframe with the count value of the genes without low expressed genes
  colnames_genes = colnames(comparison_df[1,])
  dataframe_filtered_count <- data.frame()
  for (i in 1:nrow(df_cpm_filter)) {
    core_genes=df_cpm_filter[i,1]
    coregenes=as.character(core_genes)
    ligne=subset(comparison_df,comparison_df[,1]==coregenes)
    dataframe_filtered_count <- rbind(dataframe_filtered_count,ligne)
  }

  name_filter_count <- paste0("output_filter_count_", condition1, "_vs_", condition2)
  output_filter_count_comparaison <- file.path(output_filter_count_path, paste0(name_filter_count, ".tabular")) 

  # Matrix transformation for cpm calculation
  write.table(dataframe_filtered_count, output_filter_count_comparaison, quote = FALSE, row.names = FALSE, col.names = FALSE)

}




