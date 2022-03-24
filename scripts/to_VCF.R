library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

non_VCF <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)

VCF <-non_VCF[,c("Otherinfo1", "Otherinfo2", "Otherinfo3", "Otherinfo4", "Otherinfo5", "Otherinfo6", "Otherinfo7", "Otherinfo8", "Otherinfo9")]

colnames(VCF) <- c("#CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

write.table(VCF, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


