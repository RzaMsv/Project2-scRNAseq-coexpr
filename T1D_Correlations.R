library(Seurat)
library(tidyverse)

CD4 <- Read10X(data.dir = "Negative//")
mergedNML <- CreateSeuratObject(counts = CD4, project = "T1D", min.cells = 3, min.features = 200)

# now let's get rid of low quality cells
mergedNML <- subset(mergedNML, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA < 20000)


# Extract the gene expression matrix
gene_exp <- as.matrix(GetAssayData(object = mergedNML, slot = "counts"))
rownames(gene_exp)
write.csv(gene_exp, file = "gene_expression_matrixN.csv", row.names = TRUE)
