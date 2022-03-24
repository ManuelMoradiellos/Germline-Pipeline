library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

input_table <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)
input_table_SNV <- filter(input_table, VARIANT_CLASS == "SNV")

write.table(input_table_SNV, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
