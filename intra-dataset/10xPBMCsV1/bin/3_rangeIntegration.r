args=commandArgs(T)
cutOff_enhanced<-as.numeric(args[1])
cutOff_coaccess<-as.numeric(args[2])

library(Matrix)
library(parallel)
Data<-readRDS(paste0("../output/10xpbmc5k-snap-full_enh",cutOff_enhanced,".rds"))
labels<-read.table("../input/metadata_10xpbmc_sorted.tsv",header = TRUE,sep = "\t",check.names = FALSE)
labels<-labels[labels$barcode %in% rownames(Data),]
crn<-readRDS("../input/conns_10xpbmc.rds")

positive_crn<-na.omit(crn[crn$coaccess>=cutOff_coaccess,])

func<-function(label){
    barcode<-labels[labels$label == label,]$barcode
    cell_data<-Data[rownames(Data) %in% barcode,]
    nonZeroColumnList_cell<-diff(cell_data@p)/nrow(cell_data)
    nonZeroColumnList_cell<-t(as.data.frame(nonZeroColumnList_cell))
    colnames(nonZeroColumnList_cell)<-Data@Dimnames[[2]]
    for(i in 1:length(positive_crn$coaccess)){
        peak1<-positive_crn$Peak1[i]
        peak2<-as.character(positive_crn$Peak2[i])
        peak1_split<-strsplit(peak1,split = "_")[[1]]
        peak1_left<-as.integer(as.integer(peak1_split[2])/5000)*5000+1
        peak1_right<-as.integer(as.integer(peak1_split[2])/5000)*5000+5000
        peak1<-paste0(peak1_split[1],":",peak1_left,"-",peak1_right)
#         print(peak1)
        peak2_split<-strsplit(peak2,split = "_")[[1]]
        peak2_left<-as.integer(as.integer(peak2_split[2])/5000)*5000+1
        peak2_right<-as.integer(as.integer(peak2_split[2])/5000)*5000+5000
        peak2<-paste0(peak2_split[1],":",peak2_left,"-",peak2_right)
        if(! peak1 %in% Data@Dimnames[[2]]){ 
            next
        }
        if(! peak2 %in% Data@Dimnames[[2]]){ 
                next
        }
        if(nonZeroColumnList_cell[,peak1]>=cutOff_enhanced){
            cell_data[,peak2] = 1  
        }
        if(nonZeroColumnList_cell[,peak2]>=cutOff_enhanced){
            cell_data[,peak1] = 1  
        } 
    }
    return(cell_data)
}

cl.cores <- 8
cl <- makeCluster(cl.cores-1,type = "FORK") 
results <- parLapply(cl, levels(labels$label),  func)
integrated_Data<-do.call('rbind',results)
stopCluster(cl)

integrated_Data<-integrated_Data[order(rownames(integrated_Data)),]
saveRDS(integrated_Data,file =paste0('../output/10xpbmc5k-snap-full_enh',cutOff_enhanced,"_int",cutOff_coaccess,'.rds'))

print(paste0('enh ',cutOff_enhanced,' and int ',cutOff_coaccess,' is finished!'))
