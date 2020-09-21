args=commandArgs(T)
# cutOff_enhanced<-as.numeric(args[1])
# cutOff_coaccess<-as.numeric(args[2])
cutOff_enhanced<-0.1
cutOff_coaccess<-0.4

setwd("/Users/cuizhe/project_local/2020-06-29-10xPBMC_v1_labeled/bin")
library(Matrix)
library(parallel)
library(pbapply)
Data<-readRDS(paste0("../output/10xpbmc5k-snap-full_enh",cutOff_enhanced,".rds"))
labels<-read.table("../input/pbmc_5k_atac_label_v1.txt",header = TRUE,sep = "\t",check.names = FALSE)
labels<-labels[labels$barcode %in% rownames(Data),]
crn<-readRDS("../input/conns_10xpbmc.rds")

positive_crn<-na.omit(crn[crn$coaccess>=cutOff_coaccess,])

func<-function(label){
    barcode<-labels[labels$label == label,]$barcode
    cell_data<-Data[rownames(Data) %in% barcode,]
    nonZeroColumnList_cell<-diff(cell_data@p)/nrow(cell_data)
    nonZeroColumnList_cell<-t(as.data.frame(nonZeroColumnList_cell))
    colnames(nonZeroColumnList_cell)<-Data@Dimnames[[2]]
    list1=c()
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
            list1<-c(list1,peak2)
        }
        if(nonZeroColumnList_cell[,peak2]>=cutOff_enhanced){
            list1<-c(list1,peak1)
        }

        # print(paste0(label," : ",i,"/",length(positive_crn$coaccess)))
    }
    list1<-unique(list1)
    cell_data<-as.matrix(cell_data)
    cell_data[,list1] <-1
    cell_data<-as(cell_data,"dgCMatrix")
    return(cell_data)
}

cl.cores <- 8
cl <- makeCluster(cl.cores,type = "FORK")
results <- pblapply(cl=cl, X=unique(labels$label),  FUN=func)
integrated_Data<-do.call('rbind',results)
stopCluster(cl)

integrated_Data<-integrated_Data[order(rownames(integrated_Data)),]
saveRDS(integrated_Data,file =paste0('../output/10xpbmc5k-snap-full_enh',cutOff_enhanced,"_int",cutOff_coaccess,'.rds'))

print(paste0('enh ',cutOff_enhanced,' and int ',cutOff_coaccess,' is finished!'))
