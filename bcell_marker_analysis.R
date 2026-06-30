library(tidyverse)
library(pheatmap)
library(ggrepel)
source('data_analysis.R')
source('differential_expression_analysis.R')

#focus on b-cell markers and cd20: MS4A1 (CD20), CD19, CD22, CD79A, CD79B, PAX5, BANK1, BLNK
bcell_genes <- c("MS4A1", "CD19", "CD22", "CD79A", "CD79B", "PAX5", "BANK1", "BLNK")

#get Gene_Id (ensembl id) for each symbol
bcell_ids <- my_data$Gene_Id[my_data$gene_name %in% bcell_genes]
names(bcell_ids) <- my_data$gene_name[my_data$gene_name %in% bcell_genes]
#check if they were successfully found
print(bcell_ids)

#get cpm values for each gene from step 3
bcell_cpm <- cpm_matrix[bcell_ids, ]

#relabel rows with gene symbols
rownames(bcell_cpm) <- names(bcell_ids) 
#round to one decimal place
print(round(bcell_cpm, 1))

#HEATMAP--------------------------------------------------------------
#log transformation
log_bcell_cpm <- log2(bcell_cpm + 1)

#annotation columns so pheatmap can label patient/whetheer it was pre or post
annotation_col <- data.frame(
  patient   = sample_info$patient,
  timepoint = sample_info$timepoint
)
rownames(annotation_col) <- sample_info$sample

pheatmap(log_bcell_cpm,
         scale = "row", #z-score each gene
         cluster_rows = FALSE, #keep genes in the order we listed them
         cluster_cols = FALSE, #keep pre/post pairs togethe
         annotation_col = annotation_col,
         display_numbers = round(log_bcell_cpm, 1),
         number_color = "black",
         fontsize_number = 7,
         main = "B-cell marker expression: pre vs post mosunetuzumab\n(log2 CPM, row z-scored)")

#gene expression plots with pre vs post for each patient
bcell_long <- as.data.frame(bcell_cpm) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cpm") %>%
  left_join(sample_info[, c("sample", "patient", "timepoint")], by = "sample") %>%
  mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")))

ggplot(bcell_long, aes(x = timepoint, y = cpm, group = patient, color = patient)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(title = "B-cell marker expression per patient: Pre vs Post mosunetuzumab",
       x = NULL, y = "CPM") +
  theme_minimal() +
  theme(legend.position = "bottom")

#log2 fc table with pre vs post for each patient
fc_table <- bcell_long %>%
  dplyr::select(gene, patient, timepoint, cpm) %>%
  group_by(gene, patient, timepoint) %>%
  summarise(cpm = mean(cpm), .groups = "drop") %>%
  pivot_wider(names_from = timepoint, values_from = cpm) %>%
  mutate(log2FC = log2((Post + 1) / (Pre + 1))) %>%
  dplyr::select(gene, patient, Pre, Post, log2FC)

print(fc_table)

#wide view = one row per gene, one column per patients fc
fc_wide <- fc_table %>%
  dplyr::select(gene, patient, log2FC) %>%
  pivot_wider(names_from = patient, values_from = log2FC)

print(fc_wide)


write.csv(fc_wide, "bcell_marker_log2FC.csv", row.names = FALSE)
write.csv(as.data.frame(bcell_cpm) %>% rownames_to_column("gene"),
          "bcell_marker_cpm.csv", row.names = FALSE)

#get limma results of these genes from step 3
bcell_de_results <- all_genes[all_genes$gene_name %in% bcell_genes, ]
print("limma results:")
print(bcell_de_results[, c("gene_name", "logFC", "P.Value", "adj.P.Val")])