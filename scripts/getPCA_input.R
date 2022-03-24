library(dplyr)

input_file <- snakemake@input[[1]]
output_file_2levs<- snakemake@output[["PCA2levs"]]
output_file_3levs<- snakemake@output[["PCA3levs"]]


tcga <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = NULL)

print(dim(tcga))

tcga_nondup <- tcga[!(duplicated(tcga[,c("Otherinfo3", "SAMPLE")])), c("Otherinfo3", "SAMPLE", "GT")]

print(head(tcga_nondup))

all_variants <- sort(unique(tcga_nondup$Otherinfo3))

all_samples <- sort(scan(what = "character", file = snakemake@params[["samples_ID"]]))

tcga_nondup_11 <- filter(tcga_nondup, GT == "1/1")
tcga_nondup_01 <- filter(tcga_nondup, GT == "0/1")

PCA_01 <- table(factor(tcga_nondup_01$Otherinfo3, levels = all_variants),
                          factor(tcga_nondup_01$SAMPLE, levels = all_samples)) 
PCA_11 <- table(factor(tcga_nondup_11$Otherinfo3, levels = all_variants),
                          factor(tcga_nondup_11$SAMPLE, levels = all_samples)) 

PCA_01[PCA_11 == 1] <- 2

PCA_01 <- as.data.frame.matrix(PCA_01)
PCA_01 <- PCA_01[,colSums(PCA_01) > 0]

PCA_01$variant <- row.names(PCA_01)
print(dim(PCA_01))
print(head(PCA_01))




PCA_complete <- table(factor(tcga_nondup$Otherinfo3, levels = all_variants),
                          factor(tcga_nondup$SAMPLE, levels = all_samples))
PCA_complete <- as.data.frame.matrix(PCA_complete)
PCA_complete <- PCA_complete[,colSums(PCA_complete) > 0]
PCA_complete$variant <- row.names(PCA_complete)

print(dim(PCA_complete))
print(head(PCA_complete))



write.table(PCA_01[, c("variant",grep(pattern = "^TCGA", colnames(PCA_01), value = TRUE))], file=output_file_3levs, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(PCA_complete[, c("variant",grep(pattern = "^TCGA", colnames(PCA_01), value = TRUE))], file=output_file_2levs, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(PCA_01[, c("variant",grep(pattern = "^TCGA", colnames(PCA_01), value = TRUE))], file=paste0(output_file_3levs,".match"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

