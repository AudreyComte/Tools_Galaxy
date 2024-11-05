
output_filter_count <- "/home/audrey/Downloads/count_filtered.tabular"
param5 <- "/home/audrey/Downloads/coldata2.tsv"

cts <- as.matrix(read.table(output_filter_count, header=T, row.names = 1))
coldata_read <- read.delim(param5, header=TRUE, comment.char="#", quote="")
colnames(cts) <- coldata_read[,1]

coldata <- coldata_read[,-1]
rownames(coldata) <- coldata_read[,1]
coldata$condition <- factor(coldata_read$condition)
coldata$type <- factor(coldata_read$type)

if (all(colnames(cts) %in% rownames(coldata)) & all(colnames(cts) == rownames(coldata))){
  comparison <- resultsNames(cts)
  comparison_df <- as.data.frame(comparison)[-1,]
  for (i in 1:length(comparison_df)) {
    outputname <- as.character(lapply(comparison_df[[i]], sub, pattern = "condition_", replacement = ""))
    condition <- strsplit(comparison_df[[i]], split = "_vs_")
    condition <- as.character(lapply(condition[[1]], sub, pattern = "condition_", replacement = ""))
    print(outputname)
  }
}

library(MASS)
library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(tidyr)



count_filtered_path <- "/home/audrey/Downloads/count_filtered.tabular"
coldata_path <- "/home/audrey/Documents/Bulk_RNAseq_analyses/coldata_DKVD.tsv"

# Load data
count_filtered_df <- read_delim(count_filtered_path, col_names = FALSE)

# Load metadata
coldata_df <- read_delim(coldata_path, delim ="\t")


for (col in 2:ncol(count_filtered_df)) {
  value <- count_filtered_df[1, col]
  if (value %in% coldata_path$project) {
    colnames(count_filtered_df)[col] <- value
  }
}

print(count_filtered_df)

# Charger la librairie readr
library(readr)
coldata <- "/home/audrey/Documents/Bulk_RNAseq_analyses/coldata_DKVD.tsv"

# Lire le fichier TSV dans un data frame
coldata_path <- read_delim(coldata, delim = "\t")

# Vérifier la structure de coldata_path
str(coldata_path)

count_filtered_path <- "/home/audrey/Downloads/count_filtered.tabular"
count_filtered_df <- read_delim(count_filtered_path, col_names = FALSE)
str(count_filtered_df)



# Charger la librairie readr
library(readr)

# Lire le fichier coldata_path TSV dans un data frame
coldata_path <- read_delim("/home/audrey/Documents/Bulk_RNAseq_analyses/coldata_DKVD.tsv", delim = "\t")

# Vérifier la structure de coldata_path
str(coldata_path)

# Lire le fichier count_filtered_df TSV dans un data frame
count_filtered_df <- read_delim("/home/audrey/Downloads/count_filtered.tabular", delim = " ")

# Vérifier la structure de count_filtered_df
str(count_filtered_df)

# Parcourir les colonnes de count_filtered_df à partir de la deuxième colonne
for (col in 2:ncol(count_filtered_df)) {
  # Obtenir la valeur de la première ligne de la colonne courante
  value <- count_filtered_df[1, col]
  
  # Vérifier si cette valeur existe dans la colonne 'project' de coldata_path
  if (value %in% coldata_path$project) {
    # Renommer la colonne courante avec cette valeur
    colnames(count_filtered_df)[col] <- value
  }
}

# Afficher le data frame modifié
print(count_filtered_df)



# Charger la librairie readr
library(readr)

# Lire le fichier coldata_path TSV dans un data frame
coldata_path <- read_delim("/home/audrey/Documents/Bulk_RNAseq_analyses/coldata_DKVD.tsv", delim = "\t")

# Vérifier la structure de coldata_path
str(coldata_path)

# Lire le fichier count_filtered_df TSV dans un data frame
count_filtered_df <- read_delim("/home/audrey/Downloads/count_filtered.tabular", delim = " ")

# Vérifier la structure de count_filtered_df
str(count_filtered_df)

# Parcourir les colonnes de count_filtered_df à partir de la deuxième colonne
for (col in 2:ncol(count_filtered_df)) {
  # Obtenir la valeur de la première ligne de la colonne courante
  value <- colnames(count_filtered_df)[col]
  
  # Vérifier si cette valeur existe dans la colonne 'project' de coldata_path
  if (value %in% coldata_path$project) {
    # Renommer la colonne courante avec cette valeur
    colnames(count_filtered_df)[col] <- value
  }
}

# Afficher le data frame modifié
print(count_filtered_df)










library(readr)
library(dplyr)

# Lire le DataFrame sans les noms de colonnes
count_filtered_df <- read_delim("/home/audrey/Downloads/count_filtered.tabular", col_names = FALSE)

# Copier la première ligne pour l'utiliser comme nom des colonnes
column_names <- as.character(unlist(count_filtered_df[1, 2:ncol(count_filtered_df)]))


colnames(count_filtered_df)[2:ncol(count_filtered_df)] <- column_names

# Assigner les nouveaux noms de colonnes
colnames(count_filtered_df) <- column_names

file_name <- "/home/audrey/Downloads/result.tabular"
write.table(count_filtered_df, file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Afficher le DataFrame résultant
print(count_filtered_df)
