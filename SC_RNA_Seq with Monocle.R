# Cell Clustering and DE test using monocle

# load in library and data
library(monocle)
library(HSMMSingleCell)
library(reshape2)

data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

# look at data
dim(HSMM_expr_matrix)
head(HSMM_expr_matrix)
head(HSMM_sample_sheet)
head(HSMM_gene_annotation)


# Construct CellDataSet object
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd,
                       lowerDetectionLimit=1, expressionFamily=negbinomial())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Detect genes - this is a filtering step to select genes for downstream clustering
HSMM <- detectGenes(HSMM, min_expr = 0.1)
head(fData(HSMM))

expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
length(expressed_genes)
head(pData(HSMM)) # The last column can be used to filter cells


# Pseudo-time estimation
# Selecting number of genes for ordering
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, mean_expression>10 & dispersion_empirical >= 5 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, ordering_genes)

# Dimension reduction
HSMM <- reduceDimension(HSMM, max_components=2)
# Order cells -- Pseudo-time estimation
HSMM <- orderCells(HSMM, reverse=TRUE)
plot_cell_trajectory(HSMM, color_by="Hours")

# Show a few genes
my_genes <- row.names(subset(fData(HSMM), gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="Hours")

# Select a number of marker genes for DE test
marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5", "ANPEP", "PDGFRA",
                                                        "MYOG", "TPM1", "TPM2", "MYH2", "MYH3", "NCAM1", "TNNT1", "TNNT2", "TNNC1",
                                                        "CDK1", "CDK2", "CCNB1", "CCNB2", "CCND1", "CCNA1", "ID1")))

# Differential expression test
diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr="~Media")
head(diff_test_res)

# visualize the data for a couple genes
MYOG_ID1 <- HSMM[row.names(subset(fData(HSMM),
                                  gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping="Media", ncol=2)


