library(Matrix)
buenrostro2018Data<-readRDS("../input/buenrostro2018-snap-full.rds")
labels<-read.table("../input/metadata_buenrostro2018_sorted.tsv",header = TRUE,sep = "\t",check.names = FALSE)

enhanced_buenrostro2018Data<-buenrostro2018Data
cutOff<-0.1
for(label in levels(labels$label)){
    barcode<-labels[labels$label == label,]$barcode
    cell_data<-buenrostro2018Data[rownames(buenrostro2018Data) %in% barcode,]
    nonZeroColumnList_cell<-diff(cell_data@p)/nrow(cell_data)
    candidateRange<-which(nonZeroColumnList_cell>=cutOff)
    j<-1
    for(i in candidateRange){
        enhanced_buenrostro2018Data[rownames(enhanced_buenrostro2018Data) %in% barcode,i] = 1
        print(paste0(label,":",j,"/",length(candidateRange)," finished"))
        j<-j+1
    }
}
saveRDS(enhanced_buenrostro2018Data,file = paste0('../output/buenrostro2018-snap-full_enh',cutOff,'.rds'))


