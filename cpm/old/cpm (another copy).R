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
input_file2 <- getOptionValue("--galaxy_input2")
output_file1 <- getOptionValue("--galaxy_output1")
output_file2 <- getOptionValue("--galaxy_output2")
output_file3 <- getOptionValue("--galaxy_output3")
output_file4 <- getOptionValue("--galaxy_output4")
param1 <- getOptionValue("--galaxy_param1")
param2 <- getOptionValue("--galaxy_param2")
param3 <- getOptionValue("--galaxy_param3")
param4 <- getOptionValue("--galaxy_param4")

# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)


###############################################################################################################################
library(MASS)
library(edgeR)
library(limma)

# Browse the list of file count to produce a matrix
file_list <- list(input_file1, input_file2)
file_list
output_count= output_file1
output_cpm=output_file2
output_filter_count=output_file3

# Dataframe with the gene count for selected sample
dataframe_total_count <- data.frame()

# Remove the duplicate sequenced library with error in the fastqc report

rmrun_list = as.list(strsplit(param3, ",")[[1]])
if(length(rmrun_list)!=0){
  for (i in 1:length(rmrun_list)) {
      name <- rmrun_list[[i]]
      rmrun_file <- as.numeric(grep(pattern = paste(name,"_count.txt$", sep=""), file_list))
      file_list <- file_list[-rmrun_file]
  }
}

for (i in seq_along(file_list)) {
  # Lire le fichier
  data <- read.table(file_list[i], header = FALSE, comment.char = "#", quote = "", sep = "\t")
  data

  # Copier la première colonne dans le data.frame
  if (i == 1) {
    dataframe_total_count[, 1] <- data[, 1]
    dataframe_total_count[, 2] <- data[, 2]
  } else {
    # Copier seulement la deuxième colonne dans la colonne suivante du data.frame
    dataframe_total_count <- cbind(dataframe_total_count, data[, 2])
  }
}

# write the dataframe with count value of each samples


if (length(file_list) > 0) {
  write.table(dataframe_total_count, file = output_count, sep="\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
}

#gene_name <- read.delim(file_list[[1]], header=FALSE, comment.char="#", quote="")
#dataframe_total_count <- cbind(gene_name[,1])
#for (i in 1:length(file_list)) {
  #count <- read.delim(file_list[[i]], header=FALSE, comment.char="#", quote="")
  #count[1,] <- lapply(count[1,], sub, pattern = ".*\\/.*\\/", replacement = "")
  #dataframe_total_count <- cbind(dataframe_total_count, count[,2])
#}



