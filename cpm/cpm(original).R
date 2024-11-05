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
outputfile_1 <- getOptionValue("--galaxy_output4")
outputfile_2 <- getOptionValue("--galaxy_output5")
outputfile_3 <- getOptionValue("--galaxy_output6")
outputfile_4 <- getOptionValue("--galaxy_output7")
outputfile_5 <- getOptionValue("--galaxy_output8")
outputfile_6 <- getOptionValue("--galaxy_output9")
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

#input_File2 <- input_file2
output_count <- output_file1
output_cpm <- output_file2
output_filter_count <- output_file3

output_1 <- outputfile_1
output_2 <- outputfile_2
output_3 <- outputfile_3
output_4 <- outputfile_4
output_5 <- outputfile_5
output_6 <- outputfile_6

name_file <- list(output_1, output_2, output_3, output_4, output_5, output_6)

# Browse the list of file count to produce a matrix
#file_list <- strsplit(input_file1, ", ")[[1]]

file_list <- strsplit(input_file1, ",")[[1]]

if(param3!=""){
  rmrun_list = as.list(strsplit(param3, ",")[[1]])
}
if(length(rmrun_list)!=0){
  for (i in 1:length(rmrun_list)) {
      name <- rmrun_list[[i]]
      rmrun_file <- as.numeric(grep(pattern = paste(name,"_count.txt$", sep=""), file_list))
      file_list <- file_list[-rmrun_file]
  }
}

gene_name <- read.delim(file_list[[1]], header = FALSE)

dataframe_total_count <- data.frame()

dataframe_total_count <- cbind(gene_name[,1])
for (i in 1:length(file_list)) {
  count <- read.delim(file_list[[i]], header=FALSE, comment.char="#", quote="")
  count[1,] <- lapply(count[1,], sub, pattern = ".*\\/.*\\/", replacement = "")
  dataframe_total_count <- cbind(dataframe_total_count, count[,2])
}

# write the dataframe with count value of each samples
write.table(dataframe_total_count, file = output_count, sep="\t", quote = FALSE,row.names = FALSE, col.names = FALSE)


# Matrix transformation for cpm calculation
mtx_total_count <- as.matrix(dataframe_total_count[-1,-1])
class(mtx_total_count) <- "numeric"
cpm <- cpm(mtx_total_count)
df_cpm <- as.data.frame(cpm)
df_cpm <- cbind(Gene_id=dataframe_total_count[-1,1], df_cpm)

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

df_cpm_filter <- rbind(dataframe_total_count[1,],df_cpm_filter)
write.table(df_cpm_filter, file = output_cpm, sep=" ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Dataframe with the count value of the genes without low expressed genes
colnames_genes = colnames(dataframe_total_count[1,])
dataframe_filtered_count <- data.frame()
for (i in 1:nrow(df_cpm_filter)) {
  core_genes=df_cpm_filter[i,1]
  coregenes=as.character(core_genes)
  ligne=subset(dataframe_total_count,dataframe_total_count[,1]==coregenes)
  dataframe_filtered_count <- rbind(dataframe_filtered_count,ligne)
}
# Matrix transformation for cpm calculation
write.table(dataframe_filtered_count, file = output_filter_count, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Writing on the snakemake outputfile
#cpm_filtering_output=output_file4
#writeLines(c("cpm filtering step done"), cpm_filtering_output)

data_filter_count <- read.table(output_filter_count, header = FALSE, sep = "")
#file_filter_count <- data.frame(data_filter_count[, 1], data_filter_count[, 2])
#write.table(file_filter_count, file = output_1, row.names = FALSE, col.names = FALSE, sep = " ")


# Nombre total de colonnes à combiner avec la première colonne, ajuster selon le nombre réel
number_of_colonnes <- 6  

# Boucle pour créer un fichier pour chaque paire de colonnes
for (i in 2:(number_of_colonnes + 1)) {
  # Créer un nouveau data frame avec la première colonne et la i-ème colonne
  file_filter_count <- data.frame(data_filter_count[, 1], data_filter_count[, i])
  file_filter_count <- file_filter_count[-1, ]
  # Écriture du fichier de sortie
  write.table(file_filter_count, file = name_file[[i-1]], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
}
