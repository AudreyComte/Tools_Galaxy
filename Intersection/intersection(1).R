# #################################################################
# 
#        Presence or absence of differential expressed 
#            genes in different comparisons
#              
# #################################################################

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

# Input
input_file_1 <- getOptionValue("--galaxy_input1")
input_file_2 <- getOptionValue("--galaxy_input2")
stat_input_file_1 <- getOptionValue("--galaxy_input3")
stat_input_file_2 <- getOptionValue("--galaxy_input4")
coldata <- getOptionValue("--galaxy_input5")

output_path1 <- getOptionValue("--output_path1")
output_path2 <- getOptionValue("--output_path2")

if (!dir.exists(output_path1)) {
  dir.create(output_path1, recursive = TRUE)
}

if (!dir.exists(output_path2)) {
  dir.create(output_path2, recursive = TRUE)
}

file_1 = gsub("^.*/", "", input_file_1)
file_2 = gsub("^.*/", "", input_file_2)

gene_list_1 = read.delim(file = input_file_1,
                           header = FALSE,
                           stringsAsFactors = FALSE)
                           
gene_list_2 = read.delim(file = input_file_2,
                           header = FALSE,
                           stringsAsFactors = FALSE)

## Data frame with common genes and corresponding fold change
# Remove the second column (basemean)
gene_list_1 = gene_list_1[,-2]
gene_list_2 = gene_list_2[,-2]

gene_list_1$foldchange_2=gene_list_2[,2][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1$pvalue_2=gene_list_2[,5][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1[,3]=gene_list_1[,5]
gene_list_1[,4]=gene_list_1$foldchange_2
gene_list_1[,5]=gene_list_1$pvalue_2

gene_list=gene_list_1[!is.na(gene_list_1[,4]),]
gene_list <- gene_list[-1,]
colnames(gene_list)[1] = "gene_id"
colnames(gene_list)[2] = paste("log2foldchange",file_1,sep = " ")
colnames(gene_list)[3] = paste("padj",file_1,sep = " ")
colnames(gene_list)[4] = paste("log2foldchange",file_2,sep = " ")
colnames(gene_list)[5] = paste("padj",file_2,sep = " ")

name1 <- paste0("common_genes", file_1, "_vs_", file_2)
output_file1 <- file.path(output_path1, paste0(name1, ".txt"))

write.table(gene_list[,1:5], file = output_file1,row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)

## Foldgene of the common genes in different selected comparison

stat_file_1 = gsub("^.*/", "", stat_input_file_1)
stat_file_2 = gsub("^.*/", "", stat_input_file_2)

stat_1 = read.delim(file = stat_input_file_1,
                           header = FALSE,
                           stringsAsFactors = FALSE)
                           
stat_2 = read.delim(file = stat_input_file_2,
                           header = FALSE,
                           stringsAsFactors = FALSE)

stat_1 = stat_1[,-2]
stat_2 = stat_2[,-2]

gene_list$foldchange_stat_1=stat_1[,2][match(gene_list$gene_id, stat_1[,1])]
gene_list$pvalue_1=stat_1[,5][match(gene_list$gene_id, stat_1[,1])]

gene_list$foldchange_stat_2=stat_2[,2][match(gene_list$gene_id, stat_2[,1])]
gene_list$pvalue_2=stat_2[,5][match(gene_list$gene_id, stat_2[,1])]

gene_list[,2]=gene_list$foldchange_stat_1
gene_list[,3]=gene_list$pvalue_1
gene_list[,4]=gene_list$foldchange_stat_2
gene_list[,5]=gene_list$pvalue_2

colnames(gene_list)[2] = paste("log2foldchange",stat_file_1,sep = " ")
colnames(gene_list)[3] = paste("padj",stat_file_1,sep = " ")
colnames(gene_list)[4] = paste("log2foldchange",stat_file_2,sep = " ")
colnames(gene_list)[5] = paste("padj",stat_file_2,sep = " ")

name2 <- paste0("stat_common_genes", file_1, "_vs_", file_2)
output_file2 <- file.path(output_path2, paste0(name2, ".txt"))

write.table(gene_list[,1:5], file = output_file2,row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)
