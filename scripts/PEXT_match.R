library(dplyr)

# Load parameters
input_file_str <- snakemake@input[["patho_input"]]
pext_file_input_str <- snakemake@input[["PEXT_file_input"]]
pext_file_output_str <- snakemake@output[["PEXT_file_output"]]

# Read table
patho_df <- read.table(file = input_file_str, sep = "\t", header = TRUE, row.names = NULL)
pext_df <- read.table(file = pext_file_input_str, sep = "\t", header = TRUE, row.names = NULL)




write.table(AFs_keep, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

