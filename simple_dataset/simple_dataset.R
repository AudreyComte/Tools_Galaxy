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
output_path <- getOptionValue("--output_path")


if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Écrire des documents d'exemple dans le répertoire de sortie
for (name in c("one", "two", "three")) {
  file_path <- file.path(output_path, paste0(name, ".txt"))
  writeLines(paste0(name, ", 1.3"), file_path)
}