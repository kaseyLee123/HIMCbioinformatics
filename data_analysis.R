library(tidyverse)
library(dplyr)
library(pheatmap)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma","edgeR","ggrepel"))
library(edgeR)
install.packages("ggcorrplot")
library(ggcorrplot)

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




#PCA CHART---------------------------------------

#PCA works better when data follows a normal distribution --> taking log 2 normalizes the raw CPM distribution
log_cpm_matrix <- log2(cpm_matrix + 1)

#run PCA-------------
#t() transposes it so its each row is one sample rather than the gene
#center subtracts mean of each gene across samples to remove the "baseline" so easier to find actual variance
#scale true when using raw cpm but we already scaled using log
pca <- prcomp(t(log_cpm_matrix), center = TRUE, scale. = FALSE)

#calculate variance of each PC---------
#PCA produces as much PCs as Samples but we can only plot 2 at a time, 
#so knowing the variance is important for knowing how much of the full picture the axes are showing
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

#build data frame for plotting--------
pca_df <- data.frame(
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  patient   = c("P29", "P29", "P73", "P73", "P100", "P100", "P195", "P195"),
  timepoint = c("Pre", "Post", "Pre", "Post", "Pre", "Post", "Pre", "Post")
)


#plot--------
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = patient, shape = timepoint)) +
  geom_line(aes(group = patient), color = "grey70") +
  geom_point(size = 4) +
  labs(
    x     = paste0("PC1 (", var_pct[1], "% variance)"),
    y     = paste0("PC2 (", var_pct[2], "% variance)"),
    title = "PCA of FL samples (pre vs post mosunetuzumab)"
  ) +
  theme_minimal()

print(pca_plot)

#SAMPLE CORRELATION HEATMAP---------------------------------------

cor_matrix <- cor(log_cpm_matrix, method = "spearman")

heatmap <- pheatmap(cor_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Sample-to-sample correlation (Spearman)")

print(heatmap)


