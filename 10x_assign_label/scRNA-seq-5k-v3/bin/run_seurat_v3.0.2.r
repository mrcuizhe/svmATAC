library(Seurat)
library(dplyr)
# library(reticulate)
# use_python("/Users/cuizhe/anaconda3/bin/python",required = TRUE)
setwd("/Users/cuizhe/project_local/2020-06-21-seurat-rna-seq-10x-5k-v3")
# load 5x v3 pbmc data
counts <- Read10X_h5("./tmp5k_v3.0.2/5k_pbmc_v3_filtered_feature_bc_matrix.h5")
rownames(counts) <- make.unique(rownames(counts))
rna <- CreateSeuratObject(counts = counts, assay = 'RNA', min.cells = 5, min.features = 500, project = '10x_RNA')
rna <- RenameCells(rna, add.cell.id = 'rna')
mito.features <- grep(pattern = "^MT-", x = rownames(x = rna), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = rna, slot = 'counts'))
rna$percent.mito <- percent.mito

# QC
rna <- subset(x = rna, subset = nCount_RNA > 2000 & nCount_RNA < 20000 & percent.mito < 0.2)

# preprocessing
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 100)
rna <- RunTSNE(rna, dims = 1:30)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.4, algorithm = 3)
rna <- RunUMAP(rna, graph = 'RNA_nn', metric = 'euclidean')
DimPlot(rna, reduction = "umap")

rna.markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(rna.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC))
top10 <- rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(rna, features = top10$gene) + NoLegend()

new.cluster.ids <- c(
  "Naive CD4+ T",
  'CD14+ Mono',
  'Memory CD4+',
  'NK',
  'B',
  'B',
  'NK',
  'CD8+ T',
  'FCGR3A+ Mono',
  'DC',
  'Platelet'
)

names(x = new.cluster.ids) <- levels(x = rna)
rna <- RenameIdents(object = rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

rna$celltype <- Idents(rna)
nk.cells <- subset(rna, subset = celltype == 'NK')
gzmk <- GetAssayData(nk.cells, assay = 'RNA', slot = 'data')['GZMK', ]
nk.cells$bright <- ifelse(gzmk > 1, 'NK bright', 'NK dim')
ctypes <- as.vector(rna$celltype)
names(ctypes) <- names(rna$celltype)
ctypes[Cells(nk.cells)] <- nk.cells$bright
rna <- AddMetaData(rna, metadata = ctypes, col.name = 'celltype')
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(rna, "./tmp5k_v3.0.2/pbmc_5k_v3.rds")
