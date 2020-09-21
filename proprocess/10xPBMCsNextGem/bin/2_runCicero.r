# cd scripts/
# samtools sort -t CB -@ 20 ../../input/atac_pbmc_5k_nextgem_possorted_bam.bam > ../../input/atac_pbmc_5k_nextgem_possorted_bam.sorted.bam

# python split_bam_file.py

# ./bam_preprocess.sh

library(cicero)
library(data.table)
library(Matrix)
library(proxy)
library(reshape2)
library(BuenColors)
library(umap)

path = '../tmp/count_reads_peaks_output/'
files <- list.files(path,pattern = "\\.txt$")
length(files)

#assuming tab separated values with a header    
datalist = lapply(files, function(x)fread(paste0(path,x))$V4) 
#assuming the same header/columns for all files
datafr = do.call("cbind", datalist)

df_regions = read.csv("../input/atac_pbmc_5k_nextgem_peaks.bed",
                      sep = '\t',header=FALSE,stringsAsFactors=FALSE)

peaknames = paste(df_regions$V1,df_regions$V2,df_regions$V3,sep = "_")

colnames(datafr) = sapply(strsplit(files,'\\_'),'[', 1)
rownames(datafr) = peaknames

metadata <- read.table('../input/metadata.tsv',
                         header = TRUE,
                         stringsAsFactors=FALSE,quote="",row.names=1)

mat_sparse = as(datafr, "dgTMatrix")
cicero_data = data.frame(cbind(Peak=rownames(datafr)[mat_sparse@i+1],
                           Cell=colnames(datafr)[mat_sparse@j+1],
                           Count=mat_sparse@x),stringsAsFactors = FALSE)
cicero_data$Count = as.numeric(cicero_data$Count)

input_cds <- make_atac_cds(cicero_data, binarize = TRUE)

pData(input_cds)$label = metadata[rownames(pData(input_cds)),'label']

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

set.seed(2019)
input_cds <- detectGenes(input_cds)

input_cds <- estimateSizeFactors(input_cds)

input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=15,
                        reduction_method = 'tSNE', norm_method = "none")

tsne_coords <- t(reducedDimA(input_cds))
data("human.hg19.genome")
genome_ref = human.hg19.genome
file_tss='../input/hg19/hg19-tss.bed'
row.names(tsne_coords) <- row.names(pData(input_cds))

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

conns <- run_cicero(cicero_cds, genome_ref) # Takes a few minutes to run

saveRDS(conns,file="../output/conns_10xpbmc.rds")


