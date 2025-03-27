# Send R errors to stder
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
input_file3 <- getOptionValue("--galaxy_input3")
input_file4 <- getOptionValue("--galaxy_input4")

param1 <- getOptionValue("--galaxy_param1")
param3 <- getOptionValue("--galaxy_param3")
param4 <- getOptionValue("--galaxy_param4")
param5 <- getOptionValue("--galaxy_param5")
param7 <- getOptionValue("--galaxy_param7")
param8 <- getOptionValue("--galaxy_param8")

output_names_stats <- getOptionValue("--galaxy_output_names_stats")

output_names_file_down <- getOptionValue("--galaxy_output_names_down")

output_names_file_up <- getOptionValue("--galaxy_output_names_up")



output_path1 <- getOptionValue("--output_path1")
if (!dir.exists(output_path1)) {
  dir.create(output_path1, recursive = TRUE)
}

output_path2 <- getOptionValue("--output_path2")
if (!dir.exists(output_path2)) {
  dir.create(output_path2, recursive = TRUE)
}
absolute_path2 <- normalizePath(output_path2)

output_path3 <- getOptionValue("--output_path3")
if (!dir.exists(output_path3)) {
  dir.create(output_path3, recursive = TRUE)
}
absolute_path3 <- normalizePath(output_path3)

output_path4 <- getOptionValue("--output_path4")
if (!dir.exists(output_path4)) {
  dir.create(output_path4, recursive = TRUE)
}
absolute_path4 <- normalizePath(output_path4)





# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)


# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################
#number <- 1
#file_list <- strsplit(input_file1, ",")[[1]]
#for (file_name in file_list) {
  #RDS <- file_name
  #name <- readLines(input_file2, n = number)
  #number <- number + 1
  #output_file1 <- paste0("diffexp_",name,".html")
  #rmarkdown::render(input = "/media/audrey/data/galaxy/tools/myTools/Deseq2_Report/diffexp.Rmd",
                    #output_format = "html_document",
                    #output_file  = output_file1,
                    #output_dir = output_path1,
                    #quiet = FALSE)
#}

# Lire toutes les lignes du fichier input_file2
names <- readLines(input_file2)
ref_levels <- readLines(input_file3)
mutant_levels <- readLines(input_file4)

# Liste des fichiers
file_list <- strsplit(input_file1, ",")[[1]]

file_list <- sort(file_list)

# Vérifier que la longueur de file_list et names est identique pour éviter les erreurs de correspondance
if (length(file_list) != length(names)) {
  stop("stop")
}

# Boucle pour chaque fichier dans file_list
for (i in seq_along(file_list)) {
  file_name <- file_list[i]
  RDS <- file_name
  namee <- names[i]
  ref_level <- ref_levels[i]
  mutant_level <- mutant_levels[i]

  # Génération du nom du fichier de sortie
  output_file1 <- paste0("diffexp_", namee, ".html")
  
  # Rendu du document RMarkdown
  rmarkdown::render(input = "/media/audrey/data/galaxy/tools/newTools/Deseq2_Report/diffexp.Rmd",
                    output_format = "html_document",
                    output_file  = output_file1,
                    output_dir = output_path1,
                    quiet = FALSE)
}

names_stats <- readLines(output_names_stats)
names_stats_sorted <- sort(names_stats)
writeLines(names_stats_sorted, output_names_stats)

names_file_down <- readLines(output_names_file_down)
names_file_down_sorted <- sort(names_file_down)
writeLines(names_file_down_sorted, output_names_file_down)

names_file_up <- readLines(output_names_file_up)
names_file_up_sorted <- sort(names_file_up)
writeLines(names_file_up_sorted, output_names_file_up)