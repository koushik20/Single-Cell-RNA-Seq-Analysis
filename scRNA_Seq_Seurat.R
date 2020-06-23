#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

#Loading PBMC dataset
pbmc.data <- Read10X(data.dir="C:/Users/koush/Downloads/Single-Cell-Analysis/filtered_gene_bc_matrices/mm10/")

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
# Create a Seurat object from a feature
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "10X_PBMC", assay = "RNA", min.cells = 3, min.features  = 200)
pbmc

# Examining few genes in the first thirty cells
pbmc.data[c("Xkr4", "Rp1", "Sox17"), 1:30]
# The '.' values in the matrix represent 0s (no molecules detected)

#  Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size

# Calculate the proportion of transcripts mapping to mitochondrial genes
# Calculate The Percentage Of All Counts That Belong To A Given Set Of Features(Genes)
pbmc[["percent.mito"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# visualizing QC metrics, and use these to filter cells
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well

par(mfrow = c(1, 2))
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

# We filter out cells that have unique gene counts (nFeature_RNA) over 6000 or less than
# 500 Note that > and < are used to define a'gate'.  
#-Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 10 )

# Normalizing the data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Remove unwanted sources of variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mito")

#Linear dimensional reduction - PCA
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualizing both cells and features that define the PCA
VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca")
DimPlot(object = pbmc, reduction = "pca")

# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "Cldn5")

# Scatter plot across single cells, replaces GenePlot
FeatureScatter(object = pbmc, feature1 = "Cldn5", feature2 = "PC_1")
FeatureScatter(object = pbmc, feature1 = "Cldn5", feature2 = "Ndufa4l2")

# Scatter plot across individual features, replaces CellPlot
CellScatter(object = pbmc, cell1 = "AAAGCAACAGGATCGA-1", cell2 = "ACTATCTTCGTAGATC-1")

# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("Igfbp7", "Vim", "Sox11"))
RidgePlot(object = pbmc, feature = c("Igfbp7", "Vim", "Sox11"))

# Heatmaps - exploration of heterogeneity in a dataset
DimHeatmap(object = pbmc, dim = 1, cells = 200, balanced = TRUE)

par(mfrow = c(1, 2))
DimHeatmap(object = pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# We identify 'significant' PCs as those who have a strong enrichment of low p-value features
pbmc <- JackStraw(object = pbmc, dims = 20, num.replicate = 100, verbose = FALSE)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)

JackStrawPlot(object = pbmc, dims = 1:15)

ElbowPlot(object = pbmc)

# Clustering the cells
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = TRUE)
pbmc <- FindClusters(pbmc, algorithm = 1 ,resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Non-linear dimensional reduction (tSNE)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
DimPlot(object = pbmc, reduction = "tsne")

#umap
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", split.by = "seurat_clusters")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster1.markers

VlnPlot(pbmc, features = c("Gria2", "Hmgb2"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("Pcp4", "Aldoc"), slot = "counts", log = TRUE)


FeaturePlot(pbmc, features = c("Pcp4", "Snca", "Cck", "Fabp7", "Ube2c", "Hmgb2", "Tubb5", "Foxg1"), cols = c("grey", "blue"), reduction = "tsne")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()
DoHeatmap(object = pbmc, features = top10$gene, label = TRUE)

#Demo code
library(SingleR)
align <- HumanPrimaryCellAtlasData()
align

CellsByIdentities(pbmc, cells = NULL)
pred.cell <- SingleR(test=pbmc, ref=ref, labels=ref$label, de.method="cluster")
table(pred.cell$labels)
pred.cell
