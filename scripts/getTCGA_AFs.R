library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

tcga <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)
metadata <- read.table(file = snakemake@params[["pop_info"]], sep = "\t", header = TRUE, row.names = NULL)


print(colnames(tcga))
print(colnames(metadata))

tcga$bcr_patient_barcode <- sapply(tcga$SAMPLE, function(x) substr(x,1,12))
tcga <- merge(tcga, metadata, by = "bcr_patient_barcode", all.x = TRUE)

print(tcga$bcr_patient_barcode[1:10])
print(metadata$bcr_patient_barcode[1:10])


tcga_nondup <- tcga[!(duplicated(tcga[,c("Otherinfo3", "SAMPLE")])), c("Otherinfo3", "SAMPLE", "GT", "mixed.PAM.Cluster")]

all_variants <- sort(unique(tcga_nondup$Otherinfo3))

populations <- readRDS(snakemake@params[["pop_counts"]])

get_AF <- function(df, population){
  # check population
  if (population == "ALL") {df_filt <- df
  } else {  df_filt <- filter(df, mixed.PAM.Cluster == population) }
  
  #count variants
  he_variants <- filter(df_filt, GT == "0/1")$Otherinfo3
  ho_variants <- rep(filter(df_filt, GT == "1/1")$Otherinfo3, 2)
  variants_count <- as.matrix(table(factor(c(he_variants, ho_variants), levels = all_variants)))
  
  # get AF
  if (population == "ALL") { variants_count_AF <- variants_count/(2*(sum(populations)))
  } else { variants_count_AF <- variants_count/(2*(populations[population]))}
  
  if (any(variants_count_AF[,1] > 1)) {stop()}
  return(variants_count_AF[,1])
}

#variant_AFs <- select(tcga_nondup, "Otherinfo3") %>% distinct()
#variant_AFs$AF_TCGA <- get_AF(tcga_nondup, "ALL")
#variant_AFs$AF_TCGA_EUR <- get_AF(tcga_nondup, "EUR")
#variant_AFs$AF_TCGA_ASN  <- get_AF(tcga_nondup, "ASIAN")
#variant_AFs$AF_TCGA_AFR <- get_AF(tcga_nondup, "AFR")
#variant_AFs$AF_TCGA_AMR <- get_AF(tcga_nondup, "AMR")

# get variants matrix
variant_matrix <- matrix(NA, ncol=5, nrow=length(all_variants))
variant_matrix[,1] <- get_AF(tcga_nondup, "ALL")
variant_matrix[,2] <- get_AF(tcga_nondup, "EUR")
variant_matrix[,3] <- get_AF(tcga_nondup, "ASIAN")
variant_matrix[,4] <- get_AF(tcga_nondup, "AFR")
variant_matrix[,5] <- get_AF(tcga_nondup, "AMR")

variant_df <- as.data.frame(cbind(all_variants, variant_matrix))
colnames(variant_df) <- c("Otherinfo3", "AF_TCGA","AF_TCGA_EUR", "AF_TCGA_ASN", "AF_TCGA_AFR", "AF_TCGA_AMR")

tcga <- merge(tcga, variant_df, by = "Otherinfo3", all.x = TRUE)
print(dim(tcga))


write.table(tcga, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

