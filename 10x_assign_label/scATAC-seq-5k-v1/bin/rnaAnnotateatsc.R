library(Seurat)
library(ggplot2)
set.seed(1234)
setwd("/Users/cuizhe/project_local/2020-06-27-seurat-atac-seq-10x-5k-v1")
# load atac peak counts and create gene activity matrix
peaks <- Read10X_h5("./input/atac_pbmc_5k_v1_filtered_peak_bc_matrix.h5")
activity.matrix <- CreateGeneActivityMatrix(
  peak.matrix = peaks,
  annotation.file = "./input/Homo_sapiens.GRCh37.82.gtf.gz"
)

# create object and filter cells
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = 'ATAC', project = '10x_ATAC')
pbmc.atac[['ACTIVITY']] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("./input/atac_pbmc_5k_v1_singlecell.csv", sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"

# process activity score matrix
DefaultAssay(pbmc.atac) <- 'ACTIVITY'
pbmc.atac <- FindVariableFeatures(pbmc.atac, nfeatures = 5000)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
# pbmc.atac <- RunPCA(pbmc.atac, npcs = 50)
# pbmc.atac <- FindNeighbors(pbmc.atac, reduction = 'pca', dims = 1:30)
# pbmc.atac <- RunUMAP(pbmc.atac, graph = 'RNA_nn')

# process peaks
DefaultAssay(pbmc.atac) <- "ATAC"
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
pbmc.atac <- RunLSI(pbmc.atac, n = 50, scale.max = NULL)
# pbmc.atac <- FindNeighbors(pbmc.atac, dims = 1:30, reduction = 'lsi')
# pbmc.atac <- FindClusters(pbmc.atac, graph = 'ATAC_nn', resoluton = 1)
# pbmc.atac <- RunTSNE(pbmc.atac, reduction = 'lsi', dims = 1:30, reduction.name = 'tsne.lsi', reduction.key = 'LSItSNE_')
# pbmc.atac <- RunUMAP(pbmc.atac, graph = 'ATAC_nn', reduction.name = 'umap.lsi', reduction.key = 'LSIUMAP_')
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:50)

pbmc.rna <- readRDS("./input/pbmc_5k_v3.rds")
pbmc.rna$tech <- "rna"
p1 <- DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
p1 + p2

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype, 
                                     weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

hist(pbmc.atac$prediction.score.max)
abline(v = 0.5, col = "red")

table(pbmc.atac$prediction.score.max > 0.5)

pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.5)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(pbmc.rna))  # to make the colors match
p1 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
p1 + p2

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p1 + p2

saveRDS(pbmc.atac.filtered, "./output/pbmc_5k_atac_v1.rds")
saveRDS(pbmc.atac.filtered$predicted.id,"./output/pbmc_5k_atac_label_v1.rds")
write.table(na.omit(pbmc.atac.filtered$predicted.id),"./output/pbmc_5k_atac_label_v1.txt",quote = F,col.names = F,sep="\t")