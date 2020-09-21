args=commandArgs(T)
cutOff<-as.numeric(args[1])

library(Matrix)
library(parallel)
corces2016Data<-readRDS("../input/corces2016-snap-full.rds")
labels<-read.table("../input/corces2016_barcode_metadata.tsv",header = TRUE,sep = "\t",check.names = FALSE)

labels<-labels[labels$barcode %in% rownames(corces2016Data),]

func<-function(label){
    barcode<-labels[labels$label == label,]$barcode
    cell_data<-corces2016Data[rownames(corces2016Data) %in% barcode,]
    nonZeroColumnList_cell<-diff(cell_data@p)/nrow(cell_data)
    candidateRange<-which(nonZeroColumnList_cell>=cutOff)
    for(i in candidateRange){
        cell_data[,i] = 1
    }
    return(cell_data)
}

cl.cores <- detectCores()
cl <- makeCluster(cl.cores-1,type = "FORK") 
results <- parLapply(cl, levels(labels$label),  func)
enhanced_corces2016Data<-do.call('rbind',results)
stopCluster(cl)
saveRDS(enhanced_corces2016Data,file = paste0('../output/corces2016-snap-full_enh',cutOff,'.rds'))
print(paste0(cutOff,' is finished'))
