library(dplyr)

# Load parameters
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
CPGs <- snakemake@params[["CPGs"]]
message(CPGs)

# Read table
table <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)
CPGs <- scan(file = CPGs, what = "character")

CPGs_table <- filter(table, Gene.refGene %in% CPGs)
write.table(CPGs_table, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

