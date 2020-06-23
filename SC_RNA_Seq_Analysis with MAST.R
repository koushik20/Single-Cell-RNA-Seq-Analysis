#scRNA-seq Analysis 
#Differential expression using MAST
#GSE85917 comparing the H1 & H9 embryonic stem cells

library(readxl)
sc_rna1 <- read_excel("H1-GSE85917.xlsx")
sc_rna2 <- read_excel("H9-GSE85917.xlsx")

H1 <- sc_rna1[-1]
row.names(H1) <- sc_rna1$Gene.ID
#View(H1)

H9 <- sc_rna2[-1]
row.names(H9) <- sc_rna2$Gene.ID
#View(H9)

dim(H1)
dim(H9)

head(H1)
head(H9)

colnames(H1)
colnames(H9)

library(MAST)

#combine into matrices
x0 = as.matrix(cbind(H1, H9))

#make a design vector
design = c(rep ("H1", ncol(H1)), rep("H9", ncol(H9)) )
design

#Data frame for cells & gene names to be used in MAST
cdata = data.frame(cellname=colnames(x0), celltype=design)
fdata = data.frame(genename=rownames(x0))

#computing log(TPM+1), which is required input for MAST
sizeFactor = colSums(x0)/1e6 
x1 = log(sweep(x0, 2, sizeFactor, FUN="/")+1)

#making single cell assay object
sca <- FromMatrix(x1, cdata, fdata)

#Threshold a matrix
cdr2 <- colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)
thres <- thresholdSCRNACountMatrix(assay(sca), nbins=200, min_per_bin=30)
assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
sca

#fitting Zero-inflated regression model and performing test
model = ~celltype + cngeneson
termTotest = "celltype" #Test the cell type specifically H1 & H9

#zlm is the function to do differential expression test
fit <- zlm(model, sca[1:100])
#likelihood ratio test
lrt <- lrTest(fit, termTotest)
#Results
head(lrt[,,3])

#Differential expression using SC2P
library(SC2P)

#Creating ExpressionSet objects, raw counts as inputs
x0 = round(as.matrix(cbind(H1, H9)))
design = data.frame(celltype=factor(c(rep("H1",ncol(H1)), rep("H9", ncol(H9)))))
colnames(x0) <- rownames(design)
phenoData <- new("AnnotatedDataFrame", data=design)
eset <- ExpressionSet(assayData = x0, phenoData=phenoData)
eset

#estimate two phases of gene expressions
data <- eset2Phase(eset[1:1000])

#two-phase DE test
de.sc2p <- twoPhaseDE(data, design = "celltype", test.which=1)

#Results
head(de.sc2p)

#Visualizing the top DEG's
#Top phase I DEG
geneName = topGene(de.sc2p, phase = "1")[1, "Gene.name"]
visGene(geneName, data, group.name = "celltype")

#Top phase II DEG
geneName = topGene(de.sc2p, phase = "2")[1, "Gene.name"]
visGene(geneName, data, group.name = "celltype")
