library(dplyr)

# Load parameters
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
threshold <- as.numeric(snakemake@params[["rare_threshold"]])
TCGA_threshold <- 0.01 

# Read table
AFs <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)

af_fields <- c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF",
               "AF_TCGA", "AF_TCGA_EUR", "AF_TCGA_AFR", "AF_TCGA_ASN", "AF_TCGA_AMR")
print(af_fields)

for (field in af_fields) {
  AFs[[field]] <- as.numeric(AFs[[field]])
}

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
#gnomad_fields <- grep(pattern = "^gnomAD", colnames(AFs_num), value = TRUE)
#af_fields <- c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF",
#               "AF_TCGA", "AF_TCGA_EUR", "AF_TCGA_AFR" "AF_TCGA_ASN", "AF_TCGA_AMR")
#for (field in af_fields) {
#  AFs_num[[field]] <- as.numeric(AFs_num[[field]])
#}
AFs_num_keep <- filter(AFs_num,
                      gnomAD_AF <= threshold,
                      gnomAD_AFR_AF <= threshold,
                      gnomAD_AMR_AF <= threshold,
                      gnomAD_EAS_AF <= threshold,
                      gnomAD_NFE_AF <= threshold,
                      gnomAD_SAS_AF <= threshold,
                      AF_TCGA <= TCGA_threshold,
                      AF_TCGA_EUR <= TCGA_threshold,
                      AF_TCGA_AFR <= TCGA_threshold,
                      AF_TCGA_ASN <= TCGA_threshold,
                      AF_TCGA_AMR <= TCGA_threshold)

AFs_keep <- rbind(AFs_NA_keep, AFs_num_keep)

write.table(AFs_keep, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

