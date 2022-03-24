library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
filterby <- snakemake@params[["func_refgene_vector"]]

input_table <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)
input_table_FUNC <- filter(input_table, Func.refGene %in% filterby)

write.table(input_table_FUNC, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
