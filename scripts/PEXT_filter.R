library("dplyr")
library("assertthat")

# Load files
patho_file_str <- snakemake@input[["patho_input"]]
pext_file_str <-  snakemake@input[["PEXT_file_input"]]
output_file <- snakemake@output[[1]]
PEXT_expr_threshold <- as.numeric(snakemake@params[["PEXT_expression_threshold"]])
cancer_type_equiv_str <- snakemake@params[["PEXT_tissue_equivalences"]]
TCGA_metadata_str <- snakemake@params[["TCGA_metadata"]]



# Read tables
patho_file <- read.table(patho_file_str, sep = "\t", header = TRUE)
pext_file <- read.table(pext_file_str, sep = "\t", header = TRUE)
cancer_type_equiv <- read.table(cancer_type_equiv_str, sep = ",", header = TRUE)
cancer_type_equiv_noNA <- filter(cancer_type_equiv, !(is.na(cancer_tissue)))
TCGA_metadata <- read.table(TCGA_metadata_str, sep = "\t", header = TRUE)

TCGA_metadata <- TCGA_metadata %>% rename(cancer_type = type)
TCGA_metadata[TCGA_metadata == "COAD" | TCGA_metadata == "READ"] <- "COADREAD"
cancer_types <- unique(TCGA_metadata$cancer_type)
cancer_types <- cancer_types[!(is.na(cancer_types))]

message(colnames(TCGA_metadata))
patho_file <- merge(patho_file, TCGA_metadata[,c("bcr_patient_barcode", "cancer_type")], by = "bcr_patient_barcode", all.x = TRUE)
# PEXT curation 
## Get ID.
colnames(pext_file)[1] <- "ID"
## Preprocessing annotation
pext_file <- pext_file %>% mutate(tx_annotation = substr(tx_annotation, start = 1, stop = nchar(tx_annotation)))

# Loop across variants
pext_by_variants_nested_list <- lapply(pext_file$tx_annotation, function(x) {
  ## Split by gene entry.
  a <- strsplit(x, split = "\\{|\\}")[[1]]
  # remove entries nchar == 1 (separators)
  entries <- a[!(nchar(a) == 1)]
  
  # For each genes of the variant.
  entries_list <- lapply(entries, function(entry) {
    entry_by_comma <- strsplit(entry, split = ",")[[1]] #Each measure
    splitted_list_both <- strsplit(entry_by_comma, split = ":") #Each measure (name + value)
    splitted_list <- lapply(splitted_list_both, function(splitted) return(splitted[2])) ## --> gets values
    names(splitted_list) <- lapply(splitted_list_both, function(splitted) return(splitted[1])) ## --> gets names
    return(splitted_list)
  })
  names(entries_list) <- sapply(entries_list, function(entry_list) entry_list$symbol) ## Set gene names as whole entries.
  return(entries_list)
})
names(pext_by_variants_nested_list) <- pext_file$ID #Name as variants.


# Cancer type to tissue
eqvs_list <- lapply(1:nrow(cancer_type_equiv_noNA), function(x) {
  targets <- strsplit(cancer_type_equiv_noNA[x,2], split = ";")[[1]]
  return(targets)
})
names(eqvs_list) <- cancer_type_equiv_noNA[,1]
assert_that(all(unique(c(unlist(eqvs_list))) %in% cancer_types), msg = "Cancer types not identified")


# Make cancer types list (PEXT sites per cancer type)
cancer_types_targets <- setNames(lapply(cancer_types, function(ct) {
  if (ct %in% unlist(eqvs_list)) {
    target_tissues <- names(eqvs_list)[sapply(eqvs_list, function(x) any(x == ct))]
    return(target_tissues)
  } else {
    return(NA)
  }
}), cancer_types)


# Expression file filter.
patho_file_EXP_FILT <- sapply(1:nrow(patho_file), function(x) {   # Each pathogenic variant
  ID_x <- patho_file$Otherinfo3[x]                                # Get ID
  if (!(ID_x %in% names(pext_by_variants_nested_list))) {         # If ID not in PEXT
    return(TRUE)                                                     # Keep
  } else {                                                        # Else:
    target_variant_PEXT <- pext_by_variants_nested_list[[ID_x]]      # Get PEXT variant info.
    patho_gene <- patho_file$Gene.refGene[x]                         # Get pathogenic gene.
    if (patho_gene %in% names(target_variant_PEXT)) {                # If pathogenic gene in PEXT variant entry.
      target_variant_PEXT_gene <- target_variant_PEXT[[patho_gene]]     # Get gene-PEXT information.
      if(all(is.na(cancer_types_targets[[patho_file$cancer_type[x]]]))){  # If cancer type is NA in PEXT equivalent table:
        exp_vector <- as.numeric(target_variant_PEXT_gene$mean_proportion)   # Get PEXT mean expression,
        return(ifelse(exp_vector > PEXT_expr_threshold |
                      is.na(exp_vector), TRUE, FALSE))                       # TRUE / FALSE based on expression threshold.
      } else {                                                          # Else:
        PEXT_tissues <- cancer_types_targets[[patho_file$cancer_type[x]]]  # Get PEXT tissues.
        exp_vector <- as.numeric(target_variant_PEXT_gene[PEXT_tissues])   # Expression vector.
        exp_vector <- exp_vector[!(is.na(exp_vector))]                     # Remove NAs.
        return(ifelse(mean(exp_vector) > PEXT_expr_threshold |
                      length(exp_vector) == 0, TRUE, FALSE))               # TRUE / FALSE based on expression threshold.
      }
    } else {                                                         # Else: gene not in file --> keep
      return(TRUE)
    }
  }
})


patho_file_filt <- patho_file[patho_file_EXP_FILT,]
write.table(patho_file_filt, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)



