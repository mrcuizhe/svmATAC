#Extract the header file
# cd ../input 
# ln -s /Volumes/LACIE_SHARE/data/2020-06-10-10x-scATAC-seq/atac_pbmc_5k_v1_possorted_bam/atac_pbmc_5k_v1_possorted_bam.bam ./    
# samtools view ./atac_pbmc_5k_v1_possorted_bam.bam -H > atac_v1_pbmc_5k_possorted_bam.header.sam

#Create a bam file with the barcode embedded into the read name
# cat <( cat atac_v1_pbmc_5k_possorted_bam.header.sam ) \ <( samtools view atac_pbmc_5k_v1_possorted_bam.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \ | samtools view -bS - > 10xpbmc5k.snap.bam

#sort the bam file by read name
# samtools sort -n -@ 15 -m 1G 10xpbmc5k.snap.bam -o 10xpbmc5k.snap.nsrt.bam

#Obtain snap
# sh run_snaptools.sh


#Obtain bin matrix
library(GenomicRanges)
library(SnapATAC)

metadata <- read.table('../input/metadata.tsv',
                         header = TRUE,
                         stringsAsFactors=FALSE,quote="",row.names=1)
metadata$label = as.character(metadata$label)
x.sp = createSnap(
file="../output/10xpbmc5k.snap",
sample="10xpbmc5k",
do.par = TRUE,
num.cores=10)
x.sp = x.sp[which(x.sp@barcode %in% rownames(metadata)),];

showBinSizes("../output/10xpbmc5k.snap");
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

bmat<-x.sp@bmat
colnames(bmat)<-x.sp@feature@elementMetadata@listData$name

fullRange<-read.table("../input/10xpbmc5k_5k_full.txt")

bmat_full<-matrix(ncol = length(fullRange$V1),nrow = nrow(bmat))

rownames(bmat_full)<-rownames(bmat)
colnames(bmat_full)<-fullRange$V1

for(i in 1:length(fullRange$V1)){
    range<-as.character(fullRange$V1[i])
    if(range %in% colnames(bmat)){
        bmat_full[,range]<-bmat[,range]
    } else{
        bmat_full[,range]<-0
    }
    print(paste0(i,"/",length(fullRange$V1),"is finished!"))
}

saveRDS(bmat,file="../output/10xpbmc5k-snap.rds")
saveRDS(bmat_full,file="../output/10xpbmc5k-snap-full.rds")


