library(MASS)
library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(tidyr)


# Paths for input and coldata
count_filtered_path <- "/home/audrey/Documents/count.tabular"
coldata_path <- "/home/audrey/Documents/Bulk_RNAseq_analyses/coldata.tsv"

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
  write_tsv(comparison_df, file_name)
}

print(output_file2)













library(MASS)
library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(tidyr)


# Paths for input and coldata
count_filtered_path <- "/home/audrey/Documents/count.tabular"
coldata_path <- "/home/audrey/Documents/Bulk_RNAseq_analyses/coldata.tsv"

# Load data
count_filtered_df <- read_delim(count_filtered_path, col_names = FALSE)
coldata_df <- read_delim(coldata_path)

# Update column names
colnames(count_filtered_df)[2:ncol(count_filtered_df)] <- coldata_df$project

# Detect unique conditions and create combinations
conditions <- unique(coldata_df$condition)
condition_combinations <- combn(conditions, 2, simplify = FALSE)

counter <- 0

file_list_output <- list()

# Process each condition combination
for (combo in condition_combinations) {
  counter <- counter + 1

  condition1 <- combo[1]
  condition2 <- combo[2]
  
  samples_condition1 <- coldata_df %>% filter(condition == condition1) %>% pull(project)
  samples_condition2 <- coldata_df %>% filter(condition == condition2) %>% pull(project)
  
  comparison_df <- count_filtered_df %>%
    select(X1, all_of(samples_condition1), all_of(samples_condition2)) %>%
    rename(Geneid = X1)
  
  file_name <- paste0(condition1, "_vs_", condition2, ".tabular")
  write.table(comparison_df, file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  file_list_output[[counter]] <- file_name
}

output_file2 <- file_list_output
print(output_file2)