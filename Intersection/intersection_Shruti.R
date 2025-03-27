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
input_files_up_all <- getOptionValue("--galaxy_input_files_1")
input_files_down_all <- getOptionValue("--galaxy_input_files_2")
input_files_stat_all <- getOptionValue("--galaxy_input_files_3")
input_files_names_all <- getOptionValue("--galaxy_input_files_4")

output_path1 <- getOptionValue("--output_path1")
output_path2 <- getOptionValue("--output_path2")

# Création des répertoires de sortie si non existants
if (!dir.exists(output_path1)) {
  dir.create(output_path1, recursive = TRUE)
}
if (!dir.exists(output_path2)) {
  dir.create(output_path2, recursive = TRUE)
}

# Lecture des fichiers et division des listes d'entrée
input_files_up <- strsplit(input_files_up_all, ",")[[1]]
input_files_down <- strsplit(input_files_down_all, ",")[[1]]
input_files_stat <- strsplit(input_files_stat_all, ",")[[1]]

input_files_names <- strsplit(input_files_names_all, ",")[[1]]

first_lines <- sapply(input_files_names, function(f) trimws(readLines(f, n = 1)))

print(first_lines)

add_text_file_down <- function(x) {
  paste(x, "_signif-down-regulated.txt", sep = "")
}

add_text_file_up <- function(x) {
  paste(x, "_signif-up-regulated.txt", sep = "")
}

add_text_file_stats <- function(x) {
  paste(x, "_all_gene_stats.tsv", sep = "")
}


# Appliquer la transformation à chaque élément
names_files_down <- sapply(first_lines, add_text_file_down)
names_files_up <- sapply(first_lines, add_text_file_up)
names_files_stats <- sapply(first_lines, add_text_file_stats)

# Afficher la nouvelle liste
print(names_files_stats)
print(names_files_down)
print(names_files_up)

# Combine les deux listes en une seule liste de 24 éléments
all_names <- c(names_files_down, names_files_up)
print(all_names)

combined_matrix <- as.data.frame(t(combn(all_names, 2)))
colnames(combined_matrix) <- c("File1", "File2")

# Générer les combinaisons inversées
reversed_matrix <- combined_matrix[, c("File2", "File1")]
colnames(reversed_matrix) <- c("File1", "File2")

# Fusionner les deux matrices pour inclure les inversions
full_combined_matrix <- rbind(combined_matrix, reversed_matrix)

print(full_combined_matrix)

create_verified_matrix <- function(file_names, file_list) {
  if (length(file_names) != length(file_list)) {
    stop("Error: The list of file names and the provided list have different lengths.")
  }
  matrix <- data.frame(file_names, file_list, stringsAsFactors = FALSE)
  return(matrix)
}


# Fonction pour créer une matrice vérifiée
createe_verified_matrix <- function(file_names, file_list) {
  if (length(file_names) != length(file_list)) {
    stop("Error: The list of file names and the provided list have different lengths.")
  }
  matrix <- cbind(file_names, file_list)
  return(matrix)
}

# Création des matrices pour chaque liste avec vérification des longueurs
tryCatch({
  up_matrix <- create_verified_matrix(names_files_up, input_files_up)
  print(up_matrix)
}, error = function(e) {
  print(paste("Error for UP matrix:", e$message))
})

tryCatch({
  down_matrix <- create_verified_matrix(names_files_down, input_files_down)
  print(down_matrix)
}, error = function(e) {
  print(paste("Error for DOWN matrix:", e$message))
})

tryCatch({
  stat_matrix <- create_verified_matrix(names_files_stats, input_files_stat)
  print(stat_matrix)
}, error = function(e) {
  print(paste("Error for STAT matrix:", e$message))
})


# Boucle pour chaque combinaison de fichiers
for (i in 1:nrow(full_combined_matrix)) {

  # Récupérer les noms des fichiers pour la combinaison courante
  name_input_file_1 <- full_combined_matrix[i, "File1"]
  name_input_file_2 <- full_combined_matrix[i, "File2"]
  
  # Première comparaison
  if (name_input_file_1 %in% up_matrix[, 1]) {
    index_file_1 <- which(up_matrix[, 1] == name_input_file_1)
    input_file_1 <- normalizePath(up_matrix[index_file_1, 2])
  }
  
  if (name_input_file_1 %in% down_matrix[, 1]) {
    index_file_1 <- which(down_matrix[, 1] == name_input_file_1)
    input_file_1 <- normalizePath(down_matrix[index_file_1, 2])
  }
  
  name_input_stat1 <- sub("(_signif.*|\\.txt)$", "", name_input_file_1)
  name_input_stat_file_1 <- paste0(name_input_stat1, "_all_gene_stats.tsv")
  
  if (name_input_stat_file_1 %in% stat_matrix[, 1]) {
    index_stat_file_1 <- which(stat_matrix[, 1] == name_input_stat_file_1)
    stat_input_file_1 <- normalizePath(stat_matrix[index_stat_file_1, 2])
  }
  
  # Deuxième comparaison
  if (name_input_file_2 %in% up_matrix[, 1]) {
    index_file_2 <- which(up_matrix[, 1] == name_input_file_2)
    input_file_2 <- normalizePath(up_matrix[index_file_2, 2])
  }
  
  if (name_input_file_2 %in% down_matrix[, 1]) {
    index_file_2 <- which(down_matrix[, 1] == name_input_file_2)
    input_file_2 <- normalizePath(down_matrix[index_file_2, 2])
  }
  
  name_input_stat2 <- sub("(_signif.*|\\.txt)$", "", name_input_file_2)
  name_input_stat_file_2 <- paste0(name_input_stat2, "_all_gene_stats.tsv")
  
  if (name_input_stat_file_2 %in% stat_matrix[, 1]) {
    index_stat_file_2 <- which(stat_matrix[, 1] == name_input_stat_file_2)
    stat_input_file_2 <- normalizePath(stat_matrix[index_stat_file_2, 2])
  }
  
  # Appliquer les transformations sur les données
  file_1 = gsub("^.*/", "", name_input_file_1)
  file_2 = gsub("^.*/", "", name_input_file_2)
  
  gene_list_1 = read.delim(file = input_file_1, header = FALSE, stringsAsFactors = FALSE)
  gene_list_2 = read.delim(file = input_file_2, header = FALSE, stringsAsFactors = FALSE)
  
  # Nettoyer les colonnes
  gene_list_1 = gene_list_1[, -2]
  gene_list_2 = gene_list_2[, -2]
  
  # Ajouter des colonnes pour les correspondances entre les deux listes
  gene_list_1$foldchange_2 = gene_list_2[, 2][match(gene_list_1[, 1], gene_list_2[, 1])]
  gene_list_1$pvalue_2 = gene_list_2[, 5][match(gene_list_1[, 1], gene_list_2[, 1])]
  
  gene_list_1[, 3] = gene_list_1[, 5]
  gene_list_1[, 4] = gene_list_1$foldchange_2
  gene_list_1[, 5] = gene_list_1$pvalue_2
  
  gene_list = gene_list_1[!is.na(gene_list_1[, 4]), ]
  gene_list <- gene_list[-1, ]
  
  # Nommer les colonnes
  colnames(gene_list)[1] = "gene_id"
  colnames(gene_list)[2] = paste("log2foldchange", file_1, sep = " ")
  colnames(gene_list)[3] = paste("padj", file_1, sep = " ")
  colnames(gene_list)[4] = paste("log2foldchange", file_2, sep = " ")
  colnames(gene_list)[5] = paste("padj", file_2, sep = " ")

  file_1 <- sub("\\.[^\\.]+$", "", file_1)
  file_2 <- sub("\\.[^\\.]+$", "", file_2)

  name1 <- paste0(file_1, ".txt_vs_", file_2, ".txt")
  #output_file1 <- file.path(output_path1, paste0(name1, ".txt"))

  output_file1 <- file.path(output_path1, name1)

  write.table(gene_list[,1:5], file = output_file1,row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)

  ## Foldgene of the common genes in different selected comparison

  stat_file_1 = gsub("^.*/", "", name_input_stat_file_1)
  stat_file_2 = gsub("^.*/", "", name_input_stat_file_2)

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

  name2 <- paste0("stat_common_genes_", file_1, ".txt_vs_", file_2, ".txt")
  #output_file2 <- file.path(output_path2, paste0(name2, ".txt"))

  output_file2 <- file.path(output_path2, name2)

  write.table(gene_list[,1:5], file = output_file2,row.names=FALSE,col.names = TRUE,sep="\t", quote = FALSE)

}