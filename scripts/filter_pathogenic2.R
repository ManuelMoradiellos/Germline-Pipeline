library("dplyr")

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
patho_value <- snakemake@params[["patho_value"]]
vep_file <- snakemake@params[["VEP_file"]]


#############
# Functions #
#############

vep_output_relevantcols <- function(tsv_filename){
  # Reads VEP Browser newly annotated tsv (download as .txt),
  # extracting the relevant columns to apply some filters to
  df_vep_annot <- read.table(file = tsv_filename, sep = "\t", header = TRUE)
  df_vep_annot_relevantcols <- df_vep_annot[, c('Uploaded_variation', 'Gene', 'Feature', 
                                                'LoFtool', 'SpliceAI_pred_DP_AG', 
                                                'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG',
                                                'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG',
                                                'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG',
                                                'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL')]
  colnames(df_vep_annot_relevantcols)[c(1, 2, 3)] <- c('Uploaded_variation', 'Gene.VEP', 'Feature.VEP') # Removes problematic character
  print(df_vep_annot_relevantcols$Feature.VEP[1:10])
  df_vep_annot_relevantcols$Feature.VEP <- gsub(df_vep_annot_relevantcols$Feature.VEP, pattern = '.[0-9]+$', replacement = '') # Remove version number, not useful
  print(df_vep_annot_relevantcols$Feature.VEP[1:10])
  return( as.data.frame(df_vep_annot_relevantcols %>% distinct()) ) # Features where returning duplicated rows, change type from tibble to data frame
}

last_exon_checker <- function(dataframe){
  # Obtains a subset of the dataframe keeping those rows where we have
  # NAs for the 'EXON' entry (in which exon on the gene is the variant)
  # or when that last bit its TRUE
  ## It's a bit naive as we end up with a lot of NAs, but follows Solip's directions
  dataframe$in_last_exon <- sapply(1:nrow(dataframe), function(x){
    if( dataframe[x, 'EXON'] == '-' ){
      return(NA)
    } else {
      ifelse(length(unique(strsplit(dataframe[x, 'EXON'], split = '/')[[1]])) == 1, TRUE, FALSE)
    }
  })
  return(dataframe) 
}

wholecriteria_damag_sapply <- function(dataframe, spliceai_thr = 0.8, cadd_thr = 15,
                                       loftool_thr = 0.334){
  # Applies our current criteria in the correct order: Checks if variant is 
  # located in terminal exon or if we have NA for that entry,  then checks if
  # it is splicing (if so, SpliceAI) or exonic (if so, CADD_phred or LoFtool)
  #
  # This criteria just removes variants located in the terminal exon that
  # may not truly be damaging, and keeps the rest (NAs for some predictor
  # scores, NAs and FALSE for last terminal)
  dataframe$whole_criteria_keep <- sapply(1:nrow(dataframe), function(x){
    if ( is.na(dataframe[x, 'in_last_exon']) | dataframe[x, 'in_last_exon'] == FALSE ){
      return(TRUE)
    }
    else {
      if ( dataframe[x, 'Func.refGene'] == 'splicing' ){
        if ( is.na(dataframe[x, 'spliceai_max']) ){
          return(TRUE) }  
        else {
          ifelse(dataframe[x, 'spliceai_max'] >= spliceai_thr, TRUE, FALSE) } }
      else if ( dataframe[x, 'Func.refGene'] == 'exonic' ){
        if ( !(is.na(dataframe[x, 'CADD_phred'])) ){
          ifelse( dataframe[x, 'CADD_phred'] >= cadd_thr, TRUE, FALSE) }
        else if ( !(is.na(dataframe[x, 'LoFtool'])) ){
          ifelse(dataframe[x ,'LoFtool'] < loftool_thr, TRUE, FALSE) }
        else {
          return(TRUE) } }
    }
  })
  return(dataframe)
}

important_cols_as_num <- function(dataframe){
  # Some columns are read as characters due to having '-' as NAs,
  # so we need to change their value types to perform the numeric
  # filters
  dataframe[, c("CADD_phred", "LoFtool", "SpliceAI_pred_DS_AG", 
                "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                "SpliceAI_pred_DS_DL")] <- 
    sapply(dataframe[, c("CADD_phred","LoFtool", "SpliceAI_pred_DS_AG", 
                         "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                         "SpliceAI_pred_DS_DL")], as.numeric)
  return(dataframe)
}

dataframe_vep_filterer <- function(dataframe){
  # Takes a data frame and performs various boolean  tests for # 
  # filtering, although some of them will be only for variants #
  # in the terminal exon                                       #
  
  cols_to_remove <-c('Uploaded_variation', 'Gene.VEP', 'Feature.VEP', 
                     'LoFtool', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL',
                     'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG',
                     'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL',
                     'SpliceAI_pred_SYMBOL', 'in_last_exon', 'spliceai_max',
                     'whole_criteria_keep', 'pasteID_genetrans')
  
  # 0. # Changes some column types
  dataframe <- important_cols_as_num(dataframe)
  # 1. # Checks if variants are in last exon or not
  dataframe<- last_exon_checker(dataframe)
  # 2. # Obtains the maximum SpliceAI score from each entry as #
  #      it is the one to check as recommended by the authors  #
  dataframe <- dataframe %>% 
    mutate(spliceai_max = do.call(pmax, c(select(., c("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")))))
  # 3. # Applies our filtering criteria which tests if any variant in terminal
  #      exon passes our combined criteria of SpliceAI, CADD, LoFtool
  dataframe <- wholecriteria_damag_sapply(dataframe)
  ## I don't understand why the previous function  ##
  ## returns a list and sometimes raises errors:   ##
  #dataframe$whole_criteria_keep <- unlist(dataframe$whole_criteria_keep)
  # 4. # Removes some variant types that are not of interest
  dataframe <- dataframe %>% filter(Func.refGene != 'exonic;splicing')
  # 5. # Filters out variants that do not pass our criteria
  dataframe_filt <- (dataframe %>% filter(whole_criteria_keep != FALSE))
  # 6. # Removes columns added to do the filtering as they are no further useful
  dataframe_filt <- select(dataframe_filt, -cols_to_remove)
  return(dataframe_filt)
}

#####

rare <- read.table(file = input_file, header = TRUE, sep = "\t", row.names = NULL)


print(patho_value)


# CLINVAR
if (patho_value == "clinvar") {
  write.table(rare, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# pLOF
if (patho_value == "plof") {
  
  pLOF_fields <- c("frameshift insertion", "frameshift deletion", "stopgain", "startloss", "stoploss")
  rare_pLOF <- filter(rare, (ExonicFunc.refGene %in% pLOF_fields) | (Func.refGene == "splicing"))
  
  # Work on unique variants as VarID + Gene + Transcript, create uniqueID based on that
  rare_pLOF$pasteID_genetrans <-
    paste(rare_pLOF$Otherinfo3, rare_pLOF$Gene, rare_pLOF$Feature, sep = "_")
  
  # Read VEP-Browser output keeping relevant columns and create uniqueID as before
  plof_vep_annot_nodup <- vep_output_relevantcols(vep_file)
  plof_vep_annot_nodup$pasteID_genetrans <- 
    paste(plof_vep_annot_nodup$Uploaded_variation, plof_vep_annot_nodup$Gene.VEP, plof_vep_annot_nodup$Feature.VEP, sep = "_")
  
  # Match both files and keep all rows of our input that do not match
  rareplof_matched <- merge(rare_pLOF, plof_vep_annot_nodup, by = 'pasteID_genetrans', all.x = TRUE)
  
  # Apply CADD+LoFtool/SpliceAI filters on terminal exon variants and write results
  rareplof_matched_filt <- dataframe_vep_filterer(rareplof_matched)
  write.table(rareplof_matched_filt, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Delmis
if (patho_value == "delmis") {
  write.table(rare, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
