# Functional Enrichment ======
# Load libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)
library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)
# Gene Set Enrichment ====
df <- read.csv("All_Markers.csv", header = TRUE)
original_gene_list <- df$avg_log2FC
names(original_gene_list) <- df$X
gene_list <- sort(gene_list, decreasing = TRUE)
  
gse <- gseGO(geneList = gene_list,
             ont = "BP",
             keyType = "SYMBOL",
             scoreType = "pos",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             eps = 0)

# Visualization
# Dot Plot
require(DOSE)
dotplot(gse, showCategory = 15) + ggtitle("Biological Process")

# GSEA Plot
gseaplot2(gse, geneSetID = 1, title = gse$Description[1])

# Ridge plot
ridgeplot(gse) + labs(x = "enrichment distribution")

# KEGG =========
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType ="ENTREZID",
            OrgDb = org.Hs.eg.db)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
df2 <- df[df$X %in% dedup_ids$SYMBOL,]
df2$X <- dedup_ids$ENTREZID
kegg_gene_list <- df2$avg_log2FC
names(kegg_gene_list) <- df2$X
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList = kegg_gene_list ,
               organism = kegg_organism,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType = "ncbi-geneid",
               scoreType = "pos",
               eps = 0)

# Visualization
# Dot Plot
dotplot(kk2, showCategory = 15, title = "Enriched Pathways")

# Ridge Plot 
ridgeplot(kk2) + labs(x = "enrichment distribution")

# GSEA Plot
gseaplot(kk2, title = kk2$Description[1], geneSetID = 1)

# Pathview 
melanoma <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04060", species = kegg_organism)

