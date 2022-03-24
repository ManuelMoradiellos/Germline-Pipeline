library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
threshold <- as.numeric(snakemake@params[["threshold"]])

AFs <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)

print(AFs[1:10, 140:ncol(AFs)])
AFS_filt <- filter(AFs, gnomAD_AF >= threshold, AF_TCGA >= threshold, AF_TCGA_EUR >= threshold, AF_TCGA_ASN >= threshold, AF_TCGA_AFR >= threshold, AF_TCGA_AMR >= threshold )


write.table(AFS_filt, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

