library(tidyverse)
library(ggrepel)
source('data_analysis.R')
source('differential_expression_analysis.R')
source('bcell_marker_analysis.R')

#do gene signature stuff using zscore
zscore_row <- function(x) (x - mean(x)) / sd(x)
z_bcell <- t(apply(log_bcell_cpm, 1, zscore_row))

#signature score of all 8 genes
signature_all8 <- colMeans(z_bcell)

#signature score of all 7 genes excluding cd20 to see overall b-cell vs just cd20
z_bcell_no_ms4a1 <- z_bcell[rownames(z_bcell) != "MS4A1", ]
signature_no_ms4a1 <- colMeans(z_bcell_no_ms4a1)

#cd20 alone just to compare its individual behavior
signature_ms4a1_only <- as.numeric(z_bcell["MS4A1", ])
names(signature_ms4a1_only) <- colnames(z_bcell)

sig_scores <- data.frame(
  sample                = names(signature_all8),
  signature_all8        = unname(signature_all8),
  signature_no_MS4A1    = unname(signature_no_ms4a1),
  signature_MS4A1_only  = unname(signature_ms4a1_only),
  stringsAsFactors      = FALSE
  ) %>%
  left_join(sample_info[, c("sample", "patient", "timepoint")], by = "sample") %>%
  mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")))

print(sig_scores)
write.csv(sig_scores, "bcell_signature_scores.csv", row.names = FALSE)

#plot--------------------------------------------------------------
sig_long <- sig_scores %>%
  pivot_longer(cols = c(signature_all8, signature_no_MS4A1, signature_MS4A1_only),
               names_to = "signature_type", values_to = "score") %>%
  mutate(signature_type = recode(signature_type,
                                 signature_all8        = "All 8 genes (incl. MS4A1)",
                                 signature_no_MS4A1     = "7 genes (excl. MS4A1)",
                                 signature_MS4A1_only   = "MS4A1 alone"),
         signature_type = factor(signature_type,
                                 levels = c("All 8 genes (incl. MS4A1)",
                                            "7 genes (excl. MS4A1)",
                                            "MS4A1 alone")))

print(
ggplot(sig_long, aes(x = timepoint, y = score, group = patient, color = patient)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  facet_wrap(~ signature_type) +
  labs(title = "B-cell identity signature score: Pre vs Post mosunetuzumab",
       x = NULL, y = "z-score") +
  theme_minimal() +
  theme(legend.position = "bottom")
)

#organize nicer with delta values included from pre to post
sig_combined <- sig_scores %>%
  pivot_wider(id_cols = patient,
              names_from = timepoint,
              values_from = c(signature_all8, signature_no_MS4A1, signature_MS4A1_only)) %>%
  mutate(
    delta_all8     = signature_all8_Post - signature_all8_Pre,
    delta_no_MS4A1 = signature_no_MS4A1_Post - signature_no_MS4A1_Pre,
    delta_MS4A1    = signature_MS4A1_only_Post - signature_MS4A1_only_Pre
  ) 

print(sig_combined)
write.csv(sig_combined, "bcell_signature_combined.csv", row.names = FALSE)