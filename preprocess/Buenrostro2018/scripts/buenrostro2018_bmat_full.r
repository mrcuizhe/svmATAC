library(GenomicRanges)
library(SnapATAC)

buenrostro2018_bmat<-readRDS("../input/bmat_buenrostro2018.rds")
rangeList<-readRDS("../input/rangeList_buenrostro2018.rds")
colnames(buenrostro2018_bmat)<-rangeList$name

buenrostro2018_fullRange<-read.table("../input/buenrostro2018_5k_full.txt")

buenrostro2018_bmat_full<-matrix(ncol = length(buenrostro2018_fullRange$V1),nrow = nrow(buenrostro2018_bmat))

rownames(buenrostro2018_bmat_full)<-rownames(buenrostro2018_bmat)
colnames(buenrostro2018_bmat_full)<-buenrostro2018_fullRange$V1

for(i in 1:length(buenrostro2018_fullRange$V1)){
    range<-as.character(buenrostro2018_fullRange$V1[i])
    if(range %in% colnames(buenrostro2018_bmat)){
        buenrostro2018_bmat_full[,range]<-buenrostro2018_bmat[,range]
    } else{
        buenrostro2018_bmat_full[,range]<-0
    }
    print(paste0(i,"/",length(buenrostro2018_fullRange$V1)," finished"))
}
buenrostro2018_bmat_full<-as(buenrostro2018_bmat_full,"dgCMatrix")

saveRDS(buenrostro2018_bmat_full,file="../output/buenrostro2018-snap-full.rds")





