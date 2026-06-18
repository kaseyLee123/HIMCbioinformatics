library(tidyverse)
library(dplyr)
library(DESeq2)
library(pheatmap)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma","edgeR","ggrepel"))
library(edgeR)

my_data <- read.csv("/Users/kasey/Desktop/HIMC/GSE243919_FFPE_samples_raw_read_counts.csv")
#change NA values to 0 to avoid future conflicts
my_data[is.na(my_data)] <- 0 

#create raw sample matrix
rownames(my_data) <- my_data$Gene_Id
sample_cols <- c("P29_Pre_Mosun", "P29_Post_Mosun",
                 "P73_Pre_Mosun", "P73_Post_Mosun", 
                 "P100_Pre_Mosun", "P100_Post_Mosun",
                 "P195_Pre_Mosun", "P195_Post_Mosun")
sample_matrix <- as.matrix(my_data[, sample_cols])

library_sizes <- colSums(my_data[, 2:10], na.rm = TRUE)
gene_count <- my_data %>% dplyr::count(feature == "gene")




#finding highly expressed genes by combined, pre, post
pre_samples  <- c("P29_Pre_Mosun", "P73_Pre_Mosun", "P100_Pre_Mosun", "P195_Pre_Mosun")
post_samples <- c("P29_Post_Mosun", "P73_Post_Mosun", "P100_Post_Mosun", "P195_Post_Mosun")

find_high_gene_expression <- function(dataMatrix){
  combined_totals <- rowSums(dataMatrix[,sample_cols])
  pre_totals <- rowSums(dataMatrix[,pre_samples])
  post_totals <- rowSums(dataMatrix[, post_samples])
  
  raw_comparison <- data.frame(pre_raw = pre_totals,
                               post_raw = post_totals,
                               combined_raw = combined_totals)

  results = list(top_pre_expressed = head(raw_comparison[order(raw_comparison$pre_raw, decreasing = TRUE), ], 15), 
                 top_post_expressed = head(raw_comparison[order(raw_comparison$post_raw, decreasing = TRUE), ], 15),
                 top_combined_expressed = head(raw_comparison[order(raw_comparison$combined_raw, decreasing = TRUE), ], 15))
  
  return(results)
}

#RAW---------------------------------------
raw_output <- find_high_gene_expression(sample_matrix)

top_raw_pre_expressed <- raw_output$top_pre_expressed
top_raw_post_expressed <- raw_output$top_post_expressed
top_raw_combined_expressed <- raw_output$top_combined_expressed

#CPM---------------------------------------
cpm_matrix <- cpm(sample_matrix)
cpm_output <- find_high_gene_expression(cpm_matrix)

top_cpm_pre_expressed <- cpm_output$top_pre_expressed
top_cpm_post_expressed <- cpm_output$top_post_expressed
top_cpm_combined_expressed <- cpm_output$top_combined_expressed




#create PCA chart 
#samples that cluster together should have very similar gene expression  
#while those that are far are biologically very different







