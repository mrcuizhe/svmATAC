library(GenomicRanges)
library(SnapATAC)

metadata <- read.table('../input/metadata_corces2016_sorted.tsv',
                         header = TRUE,
                         stringsAsFactors=FALSE,quote="",row.names=1)
barcode_metadata <- read.table('../tmp/corces2016_barcode_metadata.tsv',
                        header = TRUE,
                        stringsAsFactors=FALSE,quote="",row.names=1)

x.sp = createSnap(
file="../output/corces2016.snap",
sample="corces2016",
do.par = TRUE,
num.cores=10)

showBinSizes("../output/corces2016.snap");
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=10);

x.sp = makeBinary(x.sp, mat="bmat");

# system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz -O ../input/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz")

black_list = read.table('../input/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz')
black_list.gr = GRanges(black_list[,1], IRanges(black_list[,2], black_list[,3]));
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];

x.sp = filterBins(
x.sp,
low.threshold=-2,
high.threshold=2,
mat="bmat"
);

corces2016_bmat<-x.sp@bmat
colnames(corces2016_bmat)<-x.sp@feature@elementMetadata@listData$name

corces2016_fullRange<-read.table("../input/corces2016_5k_full.txt")

corces2016_bmat_full<-matrix(ncol = length(corces2016_fullRange$V1),nrow = nrow(corces2016_bmat))

rownames(corces2016_bmat_full)<-rownames(corces2016_bmat)
colnames(corces2016_bmat_full)<-corces2016_fullRange$V1

for(i in 1:length(corces2016_fullRange$V1)){
    range<-as.character(corces2016_fullRange$V1[i])
    if(range %in% colnames(corces2016_bmat)){
        corces2016_bmat_full[,range]<-corces2016_bmat[,range]
    } else{
        corces2016_bmat_full[,range]<-0
    }
}

saveRDS(corces2016_bmat,file="../output/corces2016-snap.rds")
saveRDS(corces2016_bmat_full,file="../output/corces2016-snap-full.rds")


