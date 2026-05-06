# Dimensional Reduction ================

# PCA -------------------------------------
seu <- RunPCA(seu, features = VariableFeatures(seu))
DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE) # PCA 1, most variation
DimPlot(seu, reduction = 'pca') # Threshold around 16 

ElbowPlot(seu)

# Clustering -------
seu <- FindNeighbors(seu, dims = 1:16)
seu <- FindClusters(seu, resolution = c(0.1, 0.3, 0.5))

DimPlot(seu, group.by = 'RNA_snn_res.0.3', label = TRUE)

# Metadata is original merged with seu 
nrow(metadata) = ncol(seu) #TRUE
seu$Malignant_Status <- metadata[, "Malignant_Status"]
seu$Cell_Type <- metadata[, "Cell_Type"]
seu$Tumor <- metadata[, "Tumor"]

# UMAP -------
seu <- RunUMAP(seu, dims = 1:16)
clusters <- DimPlot(seu, reduction = 'umap', label = T)

# Malignant (1 = no, 2 = yes, 0 = unresolved)
malignant <- DimPlot(seu, reduction = 'umap', group.by = 'Malignant_Status')

# Non-malignant (1=T, 2=B,3=Macro.4=Endo.,5=CAF, 6=NK)
cellt <- DimPlot(seu, reduction = 'umap', group.by = 'Cell_Type')

clusters|malignant|cell

saveRDS(seu, file = 'melanoma_seu')

