library(dplyr)

# Load parameters
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
threshold <- as.numeric(snakemake@params[["rare_threshold"]])
TCGA_threshold <- 0.01 

# Read table
AFs <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)

# Filter dataframe with bars.
AFs_NA <- filter(AFs, gnomAD_AF == "-")
AFs_NA_keep <- filter(AFs_NA,
                      AF_TCGA <= TCGA_threshold,
                      AF_TCGA_EUR <= TCGA_threshold,
                      AF_TCGA_AFR <= TCGA_threshold,
                      AF_TCGA_ASN <= TCGA_threshold,
                      AF_TCGA_AMR <= TCGA_threshold)

# Filter dataframe with numbers.
AFs_num <- filter(AFs, !(gnomAD_AF == "-"))
gnomad_fields <- grep(pattern = "^gnomAD", colnames(AFs_num), value = TRUE)
for (field in gnomad_fields) {
  AFs_num[[field]] <- as.numeric(AFs_num[[field]])
}
AFs_num_keep <- filter(AFs_num,
                      gnomAD_AF <= threshold,
                      gnomAD_AFR_AF <= threshold,
                      gnomAD_AMR_AF <= threshold,
                      gnomAD_EAS_AF <= threshold,
                      gnomAD_NFE_AF <= threshold,
                      gnomAD_SAS_AF <= threshold)

AFs_keep <- rbind(AFs_NA_keep, AFs_num_keep)


filter_AD <- sapply(AFs_keep$AD, function(x) {
    ADs_second <- as.numeric(strsplit(x, split = ",")[[1]][2])
    if (is.na(ADs_second)) {return(FALSE)}
    else {return(ADs_second >= 5)}
})

AFs_keep <- AFs_keep[filter_AD,]


write.table(AFs_keep, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

