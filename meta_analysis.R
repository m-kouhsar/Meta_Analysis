
message("Loading requiered libraries...")
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(metafor)))
set.seed(12345)

####################################################################################
# Input arguments:
#   
#   Input_Dir: A directory contains Input files 
#   
#   Input_csv_file: A csv (comma delimited) files must contains the following columns (column names are NOT case sensitive).
#             Study_name: Name of study
#             DEG_file: Full path of DEG results file for each study in the Input_Dir.
#             Gene_col: Column name for gene IDs in DEG results (default is "Gene")
#             SE_col: column name for standard error values in DEG results ("Default is "SE")
#             ES_col: column name for effect size values in DEG results (Default is "ES")
#   
#   min_gene_appearance: minimum number times each gene appears across the studies (Default is number of studies or number of rows in Input_csv_file)

# For each study, a DEG results file is required. It must be a CSV file (comma delimited) contains at least 3 columns for genes, standard errors and estimate values.

###################################################################################



args <- commandArgs(T)

Input_Dir <- trimws(args[1])
Input_csv_file <- trimws(args[2])
min_gene_appearance <- as.numeric(trimws(args[3]))
OutPrefix <- trimws(args[4])

dir.create(dirname(OutPrefix) , recursive = T)

Input_csv <- read.csv(Input_csv_file , stringsAsFactors = F , header = T)
colnames(Input_csv) <- tolower(colnames(Input_csv))

if(!all(c("study_name","deg_file","gene_col","se_col","es_col") %in% colnames(Input_csv))){
  stop("The following columns are missed in Input_csv_file:",
       setdiff(c("study_name","deg_file","gene_col","se_col","es_col") , colnames(Input_csv)))
}

Input_csv$deg_file <- paste0(Input_Dir ,"/", Input_csv$deg_file)
n_study = nrow(Input_csv)
studies <- vector(mode = "list" , length = n_study)
names(studies) <- Input_csv$study_name

col_names <- vector(mode = "numeric" , length = 3)
for(i in 1:n_study){
  
  message("Reading DEG results ",i," ",Input_csv$study_name[i]," from ",Input_csv$deg_file[i])
  
  studies[[i]] = read.csv(Input_csv$deg_file[i] , header = T ,stringsAsFactors = F)
  
  # Getting gene column
  if(tolower(Input_csv$gene_col[i]) %in% tolower(colnames(studies[[i]]))){
    
    col_names[1] = match(tolower(Input_csv$gene_col[i]) , tolower(colnames(studies[[i]])))
    
  }else if("gene" %in% tolower(colnames(studies[[i]]))){
    
    col_names[1] = match("gene" , tolower(colnames(studies[[i]])))
    
  }else{
    stop("Gene column was not found in '",Input_csv$deg_file[i],"'")
  }
  
  # Getting SE column
  if(tolower(Input_csv$se_col[i]) %in% tolower(colnames(studies[[i]]))){
    
    col_names[2] = match(tolower(Input_csv$se_col[i]) , tolower(colnames(studies[[i]])))
    
  }else if("se" %in% tolower(colnames(studies[[i]]))){
    
    col_names[2] = match("se" , tolower(colnames(studies[[i]])))
    
  }else{
    stop("SE column was not found in '",Input_csv$deg_file[i],"'")
  }
  
  # Getting ES column
  if(tolower(Input_csv$es_col[i]) %in% tolower(colnames(studies[[i]]))){
    
    col_names[3] = match(tolower(Input_csv$es_col[i]) , tolower(colnames(studies[[i]])))
    
  }else if("es" %in% tolower(colnames(studies[[i]]))){
    
    col_names[3] = match("es" , tolower(colnames(studies[[i]])))
    
  }else{
    stop("ES column was not found in '",Input_csv$deg_file[i],"'")
  }
  
  studies[[i]] <- studies[[i]][,col_names]
  colnames(studies[[i]]) <- c("Gene" , "Standard_Error" , "Effect_Size")
  rownames(studies[[i]]) <- studies[[i]][,1]
}

message("Filtering genes appear in < ",min_gene_appearance, " studies...")
# Create a list of the gene names from each study
gene_lists <- lapply(studies , function(df){df$Gene})

# Count how many times each gene appears across the studies
gene_counts <- table(unlist(gene_lists))

# Select genes that appear in at least 5 studies
common_genes <- names(gene_counts)[gene_counts >= min_gene_appearance]

studies <- lapply(studies,function(df){df[common_genes,]})
studies <- lapply(studies, na.omit)
studies <- imap(studies , function(df,name){
  colnames(df)[-1] <- paste0(colnames(df)[-1], ".", name)
  return(df)
})

data_meta.all <- Reduce(function(x,y){merge(x,y,by="Gene",all=T)},studies)


# Initialize results storage
meta_results <- data.frame(
  Gene = NA,
  Fixed_Est = NA,
  Fixed_SE = NA,
  Fixed_pval = NA,
  Random_Est = NA,
  Random_SE = NA,
  Random_pval = NA,
  tau2 = NA,       # between-study variance
  I2 = NA          # heterogeneity
)

message("Running meta-analysis...")
# Loop through each row (gene)
for (i in 1:nrow(data_meta.all)) {
  
  meta_results[i,"Gene"] = data_meta.all$Gene[i]
  try({
    # Extract values for current gene
    
    index_ES = grep("Effect_Size", colnames(data_meta.all), fixed = TRUE)
    yi <- as.numeric(data_meta.all[i,index_ES])
    
    index_SE = grep("Standard_Error", colnames(data_meta.all), fixed = TRUE)
    sei <- as.numeric(data_meta.all[i,index_SE])
    
    # Remove NAs if any
    valid <- !is.na(yi) & !is.na(sei)
    yi <- yi[valid]; sei <- sei[valid]
    
    if (length(yi) > 1) {
      # Fixed-effect model
      fe <- rma.uni(yi = yi, sei = sei, method = "FE")
      
      # Random-effects model
      re <- rma.uni(yi = yi, sei = sei, method = "REML")
      
      # Store results
      meta_results[i,-1 ] <- c(
        fe$b, fe$se, fe$pval,         # fixed
        re$b, re$se, re$pval,         # random
        re$tau2, re$I2                # heterogeneity stats
      )
    }
  })
}

index = apply(meta_results , MARGIN = 1 , function(x){all(is.na(x[-1]))})
warning("Meta analysis was failed for ",sum(index) , " genes:\n" , paste(meta_results$Gene[index] , collapse = ";"))
meta_results <- meta_results[!index , ]

meta_results$Fixed_pval_Bonf <- p.adjust(meta_results$Fixed_pval, method = "bonferroni")
meta_results$Random_pval_Bonf <- p.adjust(meta_results$Random_pval, method = "bonferroni")

meta_results$Fixed_pval_FDR <- p.adjust(meta_results$Fixed_pval, method = "fdr")
meta_results$Random_pval_FDR <- p.adjust(meta_results$Random_pval, method = "fdr")

n_fixed_sig <- sum(meta_results$Fixed_pval < 0.05, na.rm = TRUE)
n_random_sig <- sum(meta_results$Random_pval < 0.05, na.rm = TRUE)
message("Nominally-significant genes (Fixed):", n_fixed_sig)
message("Nominally-significant genes (Random):", n_random_sig)

n_fixed_sig <- sum(meta_results$Fixed_pval_Bonf < 0.05, na.rm = TRUE)
n_random_sig <- sum(meta_results$Random_pval_Bonf < 0.05, na.rm = TRUE)
message("Bonferroni-significant genes (Fixed):", n_fixed_sig)
message("Bonferroni-significant genes (Random):", n_random_sig)

n_fixed_sig <- sum(meta_results$Fixed_pval_FDR < 0.05, na.rm = TRUE)
n_random_sig <- sum(meta_results$Random_pval_FDR < 0.05, na.rm = TRUE)
message("FDR-significant genes (Fixed):", n_fixed_sig)
message("FDR-significant genes (Random):", n_random_sig)

# Combine with original matrix
data_meta.all <- data_meta.all[!index , ]
final_results <- merge.data.frame(data_meta.all, meta_results, by = "Gene")

# Save to file
write.csv(final_results, paste0(OutPrefix , ".metafor.csv"), row.names = F)

