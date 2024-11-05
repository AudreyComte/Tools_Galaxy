# How to execute this tool
# $Rscript my_r_tool.R --input input.csv --output output.csv

library(edgeR)

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
input_file1 <- getOptionValue("--input1")
output_file1 <- getOptionValue("--output1")

input_file2 <- getOptionValue("--input2")
output_file2 <- getOptionValue("--output2")


# Print options to stderr for debugging
#cat("\n input: ", input_file)
#cat("\n output: ", output_file)

# Read in the input file
inp1 <- read.csv(file= input_file1, stringsAsFactors = FALSE)
inp2 <- read.csv(file= input_file2, stringsAsFactors = FALSE)

# Changes every value in the first column to 0
inp1$V1 <- rep(0, nrow(inp1))
inp2$V1 <- rep(0, nrow(inp2))


# Write output to a new file which will be recognized by Galaxy
write.csv(inp1, file= output_file1, row.names = FALSE)
write.csv(inp2, file= output_file2, row.names = FALSE)

cat("\n success \n")
