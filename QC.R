# Set up environment =====
rm(list = ls(all.names = TRUE))
gc() # free up memory and report memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, 
        dplyr.summmarise.inform = F)

# Load libraries
library(tidyverse)
library(Seurat)


# Set input path 
path <- "/Users/Cheesy/documents/code/melanoma"
setwd(path)

set.seed(100)

# Import data -------------
df <- read.table("GSE72056_melanoma_single_cell_revised_v2.txt.gz", header = TRUE, sep = "\t")

# Making the cell column the 1st column 
unique_names <- make.unique(as.character(df[,1])) # duplicates
rownames(df) <- unique_names
df <- df[,-1]

# Split metadata from cts 
metadata <- t(df[1:3, ])
colnames(metadata) <- c("Tumor", "Malignant_Status", "Cell_Type")

# Create counts matrix 
cts <- df[4:nrow(df), ]
cts <- as.matrix(cts) 

# Create Seurat object
seu <- CreateSeuratObject(counts = cts,
                          meta.data = metadata,
                          project = "Melanoma")
seu <- SetAssayData(seu, layer = 'data', new.data = cts) # Put expression data into log-normalized data slot

#QC =======
# Feature content + Filtering 
VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2)
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
seu <- subset(seu, subset = nFeature_RNA > 1500 & nFeature_RNA < 7000 )

#Identify highly variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seu), 10)
top10_plot <- VariableFeaturePlot(seu)
LabelPoints(plot = top10_plot, points = top10, repel = TRUE)

seu <- ScaleData(seu)

