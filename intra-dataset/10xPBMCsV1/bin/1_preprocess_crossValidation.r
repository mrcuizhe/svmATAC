setwd("/Users/cuizhe/Documents/GitHub/svmATAC/intra-dataset/10xPBMCsV1/bin")
label_10xpbmc<-read.csv('../input/metadata.tsv',sep = "\t",header = TRUE)

data<-readRDS("../input/10xpbmc5k-snap-full.rds")

label_10xpbmc<-label_10xpbmc[label_10xpbmc$barcode %in% rownames(data),]

label_10xpbmc<-label_10xpbmc[order(label_10xpbmc[,1]),]

write.table(label_10xpbmc,file="../input/metadata_10xpbmc_sorted.tsv",quote=FALSE,sep="\t",row.names = FALSE)

Cross_Validation <- function(LabelsPath, col_Index = 1, OutputDir){
  "
  Cross_Validation
  Function returns train and test indices for 5 folds stratified across unique cell populations,
  also filter out cell populations with less than 10 cells.
  It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.
  Parameters
  ----------
  LabelsPath : Cell population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "

  Labels <- as.matrix(LabelsPath)
  Labels <- as.vector(Labels[,col_Index])

  Removed_classes <- !(table(Labels) > 10)
  Cells_to_Keep <- !(is.element(Labels,names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep]

  # Getting training and testing Folds
  library(rBayesianOptimization)
  n_folds = 5
  Folds <- KFold(Labels,nfolds = n_folds, stratified = TRUE)
  Test_Folds <- c(n_folds:1)
  Train_Idx <- list()
  Test_Idx <- list()
  for (i in c(1:length(Folds))){
    Temp_Folds <- Folds
    Temp_Folds[Test_Folds[i]] <- NULL
#     print(Temp_Folds[Test_Folds[i]])
    Train_Idx[i] <- list(unlist(Temp_Folds))
#     print(Train_Idx[i])
    Test_Idx[i] <- Folds[Test_Folds[i]]
#     print(Test_Idx[i])
  }
  remove(Temp_Folds,i,Folds)
#   print(Train_Idx)
#   print(Test_Idx)
#   print(col_Index)
  save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,file = paste0(OutputDir, '/CV_folds.RData'))
}

Cross_Validation(label_10xpbmc, 2, "../tmp")


