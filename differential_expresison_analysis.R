library(tidyverse)
library(edgeR)
library(limma)
library(ggrepel)

source('data_analysis.R')

sample_info <- data.frame(
  sample    = sample_cols,
  patient   = c("P29", "P29", "P73", "P73", "P100", "P100", "P195", "P195"),
  timepoint = c("Pre", "Post", "Pre", "Post", "Pre", "Post", "Pre", "Post")
)
#4th column combines patient+timepoint for plot labeling later
sample_info$pat_time <- paste(sample_info$patient, sample_info$timepoint, sep = "_")

DGE <- DGEList(sample_matrix, samples = sample_info$sample, group = sample_info$timepoint)

#identifies genes that have enough counts (true/false)
keep.exprs <- filterByExpr(DGE)
#filters out the lowly expressed ones and recalculates library size
DGE <- DGE[keep.exprs, keep.lib.sizes = FALSE]

#TMM normalization
DGE <- calcNormFactors(DGE, method = "TMM")

#MDS PLOT---------------------------------------------
mds <- plotMDS(DGE, top = 500, plot = FALSE, gene.selection = "common")
mds_dataframe <- cbind(sample_info$sample, sample_info$pat_time, mds$x, mds$y)
colnames(mds_dataframe) <- c("sample", "pat_time", "pc1", "pc2")
mds_dataframe <- as_tibble(mds_dataframe) %>% mutate(across(c(pc1, pc2), as.numeric))
mds_plot <- mds_dataframe %>% ggplot(aes(x = pc1, y = pc2, color = pat_time)) +
  geom_point(size = 3) +
  xlab("Principal coordinate 1") +
  ylab("Principal coordinate 2") +
  geom_text_repel(aes(label = sample), size = 3) +
  ggtitle("MDS plot: pre vs post mosunetuzumab")
print(mds_plot)

#creates 2-factor(patient + timepoint) design matrix for more complex stuff 
patient   <- factor(sample_info$patient)
timepoint <- factor(sample_info$timepoint, levels = c("Pre", "Post"))
design    <- model.matrix(~ patient + timepoint)

#assigns quality weights to samples so outliers get weighed down
vwts <- voomWithQualityWeights(DGE, design = design, normalize.method = "none", plot = TRUE)
fit <- lmFit(vwts, design)
fit <- eBayes(fit, robust = TRUE)

#prints out a a
print(summary(decideTests(fit, adjust.method = "fdr", p.value = 0.05)))
volcanoplot(fit, coef = "timepointPost", highlight = 10)

#gets 10 most significant DE genes sorted by p-values  
top10 <- topTable(fit, adjust = "BH", coef = "timepointPost", resort.by = "P")
print(top10)

#full results table
all_genes <- topTable(fit, adjust = "BH", coef = "timepointPost",
                      p.value = 1, number = Inf, resort.by = "P")

#shows distribution of p-values 
print(summary(all_genes$adj.P.Val))
#how big expression differences are from pre and post
print(summary(all_genes$logFC))

#VOLCANO PLOT--------------------------------------------
all_genes$gene_id <- rownames(all_genes)
all_genes$sig <- ifelse(all_genes$adj.P.Val < 0.2 & all_genes$logFC > 0.5,  "Up",
                        ifelse(all_genes$adj.P.Val < 0.2 & all_genes$logFC < -0.5, "Down", "NS"))

ggplot(all_genes, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_text_repel(data = subset(all_genes, adj.P.Val < 0.2 & abs(logFC) > 1),
                  aes(label = gene_id), size = 2.5, max.overlaps = 20) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot: post vs pre mosunetuzumab",
       x = "Log2 fold change",
       y = "-log10(p-value)",
       color = "Expression") +
  theme_minimal()

#ranked gene list for most upregulated to most downregulated by fold change
ranked_genes <- all_genes[order(all_genes$logFC, decreasing = TRUE), ]
ranked_genes$rank <- 1:nrow(ranked_genes)
print(head(ranked_genes, 20))
print(tail(ranked_genes, 20))

#significantly upregulated genes (more relaxed thresholds)
sig_up <- all_genes[all_genes$adj.P.Val < 0.2 & all_genes$logFC > 0.5, ]
sig_up <- sig_up[order(sig_up$logFC, decreasing = TRUE), ]
print(paste("Number of upregulated genes:", nrow(sig_up)))
print(sig_up)

#significantly downregulated genes (more relaxed thresholds)
sig_down <- all_genes[all_genes$adj.P.Val < 0.2 & all_genes$logFC < -0.5, ]
sig_down <- sig_down[order(sig_down$logFC), ]
print(paste("Number of downregulated genes:", nrow(sig_down)))
print(sig_down)
