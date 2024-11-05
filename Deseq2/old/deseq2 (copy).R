# #################################################################
# 
#                               DESeq2
#                          
# #################################################################


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
input <- getOptionValue("--galaxy_input")

output_1 <- getOptionValue("--galaxy_output1")
output_2 <- getOptionValue("--galaxy_output2")
output_2 <- getOptionValue("--galaxy_output2")

thread_param <- getOptionValue("--galaxy_param_thread")
param1_ref_level <- getOptionValue("--galaxy_param1")
param2_coldata <- getOptionValue("--galaxy_param2")
#param3_rmproj_list <- getOptionValue("--galaxy_param3")

# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)


########################################################################################################################################################################################################################################################################################################################################################################################################


library(DESeq2)
library(readr)

parallel <- FALSE
if (thread_param > 1) {
    library("BiocParallel")
    setup parallelization
    register(MulticoreParam(thread_param))
    parallel <- TRUE
}
# Loading the parameters
# nom d'une condition : "dorsolized_torso"
ref_level = param1_ref_level
normalized_counts_file = output_1

# Rename column name of the count matrix as coldata
# colData and countData must have the same sample order
cts <- as.matrix(read.table(input, header=T, row.names = 1))
coldata_read <- read.delim(param2_coldata, header=TRUE, comment.char="#", quote="")
colnames(cts) <- coldata_read[,1]

coldata <- coldata_read[,-1]
rownames(coldata) <- coldata_read[,1]
coldata$condition <- factor(coldata_read$condition)
coldata$type <- factor(coldata_read$type)

#rmproj_list = as.list(strsplit(param3_rmproj_list, ",")[[1]])


#if(length(rmproj_list)!=0){
  #for (i in 1:length(rmproj_list)) {
      #name <- rmproj_list[[i]]
      #coldata <- coldata[-match((name), table = rownames(coldata)), ]
  #}
#}

# Check that sample names match in both files
if (all(colnames(cts) %in% rownames(coldata)) & all(colnames(cts) == rownames(coldata))){
  # Create the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
} else {
  print("sample names doesn't match in both files")
}

# Remove uninformative columns (to do when filter not already done with the CPM threshold)
#dds <- dds[ rowSums(counts(dds)) > 10, ]

# Specifying the reference level
dds$condition <- relevel(dds$condition, ref = ref_level)

# DESeq : Normalization and preprocessing (counts divided by sample-specific size factors
# determined by median ratio of gene counts relative to geometric mean per gene)
dds <- DESeq(dds, parallel = parallel)
# To save the object in a file for later use
all_rds = output_2
saveRDS(dds, file = all_rds)

# Already done in the DESeq function
dds <- estimateSizeFactors( dds)
print(sizeFactors(dds))
# Save the normalized data matrix
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file = normalized_counts_file, sep="\t", quote=F, col.names=NA)

# File to complete the snakemake rule
deseq2_init_output = output_3
writeLines(c("deseq2_init step done"), deseq2_init_output)
