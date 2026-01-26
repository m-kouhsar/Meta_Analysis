
message("Loading requiered libraries...")
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(metafor)))
suppressWarnings(suppressMessages(library(ggplot2)))
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

data_meta.all <- data.frame()

for (study_name in names(studies)) {
  
  # Extract the dataframe for this study
  temp_df <- studies[[study_name]]
  
  # Add a column identifying the study
  temp_df$study <- study_name
  
  # Bind it to the master dataframe
  data_meta.all <- rbind(data_meta.all, temp_df)
}


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
for (i in 1:length(common_genes)) {
  
  #print(i)
  dat <- subset(data_meta.all, Gene == common_genes[i])
  meta_results[i,"Gene"] = common_genes[i]
  try({
    
    if ((sum(!is.na(dat$Standard_Error) )> 1) & (sum(!is.na(dat$Effect_Size) )> 1)) {
      # Fixed-effect model
      fe <- rma.uni(yi = Effect_Size, sei = Standard_Error, method = "FE" , data = dat , slab = study)
      
      # Random-effects model
      re <- rma.uni(yi = Effect_Size, sei = Standard_Error, method = "REML", data = dat , slab = study)
      
      # Store results
      meta_results[i,-1 ] <- c(
        fe$b, fe$se, fe$pval,         # fixed
        re$b, re$se, re$pval,         # random
        re$tau2, re$I2                # heterogeneity stats
      )
    }else{
      message("Warning! Not enaough data for gene ",i,": ",common_genes[i])
    }
    
  })
}

index = apply(meta_results , MARGIN = 1 , function(x){all(is.na(x[-1]))})
if(sum(index)>0)
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

# Save to file
write.csv(meta_results, paste0(OutPrefix , ".metafor.csv"), row.names = F)

############################################################################################
#             Visualization
############################################################################################
PlotData <- data.frame(Gene = meta_results$Gene, Random_Pval = meta_results$Random_pval , Fixed_Pval = meta_results$Fixed_pval)

pdf(file = paste0(OutPrefix , ".metafor.plots.pdf") , width = 10,height = 8)
# 3. Create the plot
ggplot(PlotData, aes(x = Random_Pval)) + 
  # Histogram: Note the y = after_stat(density) argument
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightblue", 
                 color = "black") +
  # Density line: Overlay the line
  geom_density(color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Histogram of P-value in Random effect model", x = "P-value", y = "Density")

ggplot(PlotData, aes(x = Fixed_Pval)) + 
  # Histogram: Note the y = after_stat(density) argument
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightblue", 
                 color = "black") +
  # Density line: Overlay the line
  geom_density(color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Histogram of P-value in Fixed effect model", x = "P-value", y = "Density")

meta_results <- meta_results[order(meta_results$Random_pval , decreasing = F),]
top_genes <- meta_results$Gene[1:10]


par(mfrow = c(3, 2), mar = c(4, 4, 2, 2)) 

# 2. Loop through the top genes
for (g in top_genes) {
  
  # A. Subset data for this specific gene
  dat <- subset(data_meta.all, Gene == g)
  
  # B. Run the Random Effects Model
  # Note: We use 'sei' because you have Standard Error (seTE)
  res_re <- rma(yi = Effect_Size, sei = Standard_Error, data = dat, slab = study)
  
  # C. Run Fixed/Common Effect Model (optional, to match your image)
  res_fe <- rma(yi = Effect_Size, sei = Standard_Error, data = dat, method = "FE")
  
  k <- nrow(dat)
  # CHANGE 2: Manually set ylim
  # c(bottom_limit, top_limit)
  # -2.5: Extends the bottom to make room for the extra diamond at row -2
  # k + 3: Keeps the standard top spacing for headers
  forest(res_re,
         header = TRUE,
         main = g,
         xlab = "Effect Size",
         cex = 0.95,
         mlab = "Random Effects",
         ylim = c(-2.5, k + 3) 
  )
  
  # E. Add the Common Effect Diamond (optional)
  addpoly(res_fe, row = -2, mlab = "Common Effect")
}

# Reset layout
par(mfrow = c(1, 1))

graphics.off()








