library(dplyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
patho_value <- snakemake@params[["patho_value"]]
delmis_method <- snakemake@params[["delmis_method"]]
print(delmis_method)
rare <- read.table(file = input_file, header = TRUE, sep = "\t", row.names = NULL)


print(patho_value)


# CLINVAR
if (patho_value == "clinvar") {
  pathogenic <- c("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic", "risk_factor", "_risk_factor")
  benign <- c("uncertain_significance", "conflicting_interpretations_of_pathogenicity", "likely_benign", "benign", "benign/likely_benign")
  rare$is_patho <- sapply(1:nrow(rare), function(x) {
				  clinvar_tags <- strsplit(rare$CLIN_SIG[x], split = ",")[[1]]
				  return(any(clinvar_tags %in% pathogenic))})
  rare_patho_clinvar <- filter(rare, is_patho == TRUE) 

  write.table(rare_patho_clinvar[,1:(ncol(rare_patho_clinvar)-1)], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# pLOF
if (patho_value == "plof") {
  pLOF_fields <- c("frameshift insertion", "frameshift deletion", "stopgain", "startloss", "stoploss")
  rare_pLOF <- filter(rare, (ExonicFunc.refGene %in% pLOF_fields) | (Func.refGene == "splicing"))
  write.table(rare_pLOF, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Delmis
if (patho_value == "delmis") {
  rare <- filter(rare, ExonicFunc.refGene == "nonsynonymous SNV")
  if (delmis_method == "MetaLR") {
    rare_delmis <- filter(rare, MetaLR_pred == "D")
  } else if (delmis_method == "all") { 
    pass_matrix <- matrix(nrow = nrow(rare), ncol = 5)
    rare$CADD_phred <- as.numeric(rare$CADD_phred)
    pass_matrix[,1] <- rare$MetaLR_pred == "D"
    pass_matrix[,2] <- rare$SIFT_pred == "D"
    pass_matrix[,3] <- rare$Polyphen2_HVAR_pred == "D"
    pass_matrix[,4] <- rare$MutationTaster_pred == "D"
    pass_matrix[,5] <- sapply(rare$CADD_phred, function(x) {
      numX <- as.numeric(x)
      if (is.na(numX)) {return(FALSE)}
      else {return(numX >= 20)}
    })
				    
    indexing_vector <- rowSums(pass_matrix) == 5
    rare_delmis <- rare[indexing_vector,]
  }
  write.table(rare_delmis, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
