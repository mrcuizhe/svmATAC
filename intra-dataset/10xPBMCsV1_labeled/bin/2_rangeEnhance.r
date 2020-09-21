
args=commandArgs(T)
cutOff<-as.numeric(args[1])

setwd("/home/gaoy1/data/cuizhe/project/2020-06-29-10xPBMC_v1_labeled/bin")
library(Matrix)
library(parallel)
Data<-readRDS("../input/10xpbmc5k-snap-full.rds")
labels<-read.table("../input/pbmc_5k_atac_label_v1.txt",header = TRUE,sep = "\t",check.names = FALSE)
Data<-as(Data, "dgCMatrix")

labels<-labels[labels$barcode %in% rownames(Data),]

func<-function(label){
    barcode<-labels[labels$label == label,]$barcode
    cell_data<-Data[rownames(Data) %in% barcode,]
    nonZeroColumnList_cell<-diff(cell_data@p)/nrow(cell_data)
    candidateRange<-which(nonZeroColumnList_cell>=cutOff)
    for(i in candidateRange){
        cell_data[,i] = 1
    }
    return(cell_data)
}


cl.cores<-16
cl <- makeCluster(cl.cores,type = "FORK") 
results <- parLapply(cl, unique(labels$label),  func)
enhanced_Data<-do.call('rbind',results)
stopCluster(cl)
saveRDS(enhanced_Data,file = paste0('../output/10xpbmc5k-snap-full_enh',cutOff,'.rds'))
print(paste0(cutOff,' is finished'))
