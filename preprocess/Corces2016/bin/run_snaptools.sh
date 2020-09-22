#!/bin/bash

# samtools merge -@ 15 ../tmp/corces2016.nondup.merged.bam -b ../tmp/list_bamfiles.txt
# samtools sort -@ 15 -n ../tmp/corces2016.nondup.merged.bam -o ../tmp/corces2016.nondup.merged.sorted.bam

snaptools snap-pre  \
    --input-file=../tmp/corces2016.nondup.merged.sorted.bam  \
    --output-snap=../output/corces2016.snap  \
    --genome-name=hg19  \
    --genome-size=../input/hg19/hg19.chrom.sizes  \
    --min-mapq=30  \
    --min-flen=0  \
    --max-flen=1000  \
    --keep-chrm=TRUE  \
    --keep-single=False  \
    --keep-secondary=False  \
    --overwrite=True  \
    --max-num=1000000  \
    --min-cov=100  \
    --verbose=True

snaptools snap-add-bmat    \
    --snap-file=../output/corces2016.snap \
    --bin-size-list 5000    \
    --verbose=True

