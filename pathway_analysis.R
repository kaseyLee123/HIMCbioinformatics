library(tidyverse)
library(ggrepel)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "fgsea"))
install.packages("gprofiler2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(gprofiler2)

if (!require("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
library(ReactomePA)

source('data_analysis.R')
source('differential_expression_analysis.R')

#convert gene symbols to entrez ids (required for cluster profiling)
gene_list <- bitr(all_genes$gene_name,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)

all_genes_entrez <- merge(all_genes, gene_list,
                          by.x = "gene_name", by.y = "SYMBOL",
                          all.x = FALSE)

#GO vs KEGG vs REACTOME
#GO:good for initial overview, classifies genes by role, location, process
#KEGG: highlights signaling cascades and metabolic routes
#REACTOME: very human-centric, good for detailing molecular reactions

#OVER REPRESENTATION ANAYLYSIS--------------------------------------
#sees if any of the genes are in a specific pathway more than you would expect/more than chance

#takes just the entrezid of the significant genes (according to p val and logfc)
sig_genes <- all_genes_entrez$ENTREZID[all_genes_entrez$adj.P.Val < 0.2 & 
                                         abs(all_genes_entrez$logFC) > 0.5]
background_genes <- all_genes_entrez$ENTREZID

#go
go_ora <- enrichGO(gene          = sig_genes, 
                   universe      = background_genes,
                   OrgDb         = org.Hs.eg.db, #human gene database
                   ont           = "BP", #test biological process pathways specifically
                   pAdjustMethod = "BH", #fdr correction
                   pvalueCutoff  = 0.05, 
                   qvalueCutoff  = 0.2,
                   readable      = TRUE) #show gene symbols not entrezids

print(summary(go_ora))
print(dotplot(go_ora, showCategory = 20, title = "GO Biological Process enrichment"))

#kegg
kegg_ora <- enrichKEGG(gene          = sig_genes,
                       universe      = background_genes,
                       organism      = "hsa", #human
                       pAdjustMethod = "BH", #fdr correction
                       pvalueCutoff  = 0.05)

print(summary(kegg_ora))
print(dotplot(kegg_ora, showCategory = 20, title = "KEGG pathway enrichment"))

#tests reactome pathways
reactome_ora <- enrichPathway(gene          = sig_genes,
                              universe      = background_genes,
                              organism      = "human",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable      = TRUE)
print(summary(reactome_ora))
print(dotplot(reactome_ora, showCategory = 20, title = "Reactome pathway enrichment"))

#GENE SET ENRICHMENT ANALYSIS----------------------------------------------------
#shows whether genes in a specific pathway tend to cluster at the top or bottom together or are randomly dispersed

ranked_list <- all_genes_entrez$logFC
names(ranked_list) <- all_genes_entrez$ENTREZID
#sorted from upregulated to downregulated
ranked_list <- sort(ranked_list, decreasing = TRUE)
#removes duplicates just in case
ranked_list <- ranked_list[!duplicated(names(ranked_list))]

#go
gsea_go <- gseGO(geneList     = ranked_list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP", #biological process
                 minGSSize    = 10, #minimum 10 genes for pathway to be tested
                 maxGSSize    = 500, #maximum 500 genes for pathway to be tested
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

print(summary(gsea_go))
print(dotplot(gsea_go, showCategory = 20, split = ".sign",
              title = "GSEA GO Biological Process") +
        facet_grid(. ~ .sign))

#prints plot of most significant go pathway according to default (p-value)
if (nrow(as.data.frame(gsea_go)) > 0) {
  print(gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1]))
}

#kegg
gsea_kegg <- gseKEGG(geneList     = ranked_list,
                     organism     = "hsa",
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)

print(summary(gsea_kegg))
print(dotplot(gsea_kegg, showCategory = 20, split = ".sign",
              title = "GSEA KEGG pathways") +
        facet_grid(. ~ .sign))

#reactome
gsea_reactome <- gsePathway(geneList     = ranked_list,
                            organism     = "human",
                            minGSSize    = 10,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = FALSE)
print(summary(gsea_reactome))
print(dotplot(gsea_reactome, showCategory = 20, split = ".sign", title = "GSEA Reactome pathways") +
        facet_grid(. ~ .sign))

#gprofiler 2 automatically tests against GO KEGG and Reactome databases all at once
gp_results <- gost(query             = all_genes$gene_name[all_genes$adj.P.Val < 0.2 &
                                                             abs(all_genes$logFC) > 0.5],
                   organism          = "hsapiens",
                   ordered_query     = FALSE,
                   correction_method = "fdr",
                   sources           = c("GO:BP", "KEGG", "REAC"))

#results in a manhattan plots
print(gostplot(gp_results, capped = TRUE, interactive = FALSE))

gp_table <- gp_results$result[order(gp_results$result$p_value), ]
print(head(gp_table[, c("source", "term_name", "p_value", "intersection_size")], 30))


write.csv(as.data.frame(go_ora),   "go_ora_results.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_ora), "kegg_ora_results.csv", row.names = FALSE)


# fix list columns before saving gprofiler results
gp_table_clean <- gp_table

# convert any list columns to text automatically
for (col in colnames(gp_table_clean)) {
  if (is.list(gp_table_clean[[col]])) {
    gp_table_clean[[col]] <- sapply(gp_table_clean[[col]], paste, collapse = ",")
  }
}

write.csv(gp_table_clean, "gprofiler_results.csv", row.names = FALSE)
write.csv(gp_table_clean, "gprofiler_results.csv", row.names = FALSE)