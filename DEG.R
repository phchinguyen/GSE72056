library(dplyr)

# DEGs ========
mlm_seu <- readRDS('melanoma_seu')

# DefaultAssy must be RNA
DefaultAssay(mlm_seu) 

# FindAllMarkers =====
all_markers <- FindAllMarkers(mlm_seu, logfc.threshold = 0.5,
               min.pct = 0.25,
               only.pos = TRUE)
top5_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%     
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)  
# Able to verify the clustering correlates to the cell types/ subpopulation divided by the authpr

all_markers <- all_markers[, 1:5]

#FindMarkers ====
# Malignant vs. Non-Malignant
FindMarkers()
FeaturePlot(seu, features = c(""), min.cutoff = 'q10')

# Prepare file for enrichment analysis 
write.csv(all_markers, "All_Markers.csv", row.names = TRUE)
