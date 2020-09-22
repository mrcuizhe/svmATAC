#!/bin/bash

bampath=../../tmp/10xpbmc_bam_sorted_nondup
# dirlist=(`ls $bampath*.bam`)
# echo ${dirlist[$LSB_JOBINDEX-1]}
mkdir -p ../../tmp/count_reads_peaks_output
for file in $bampath/*.bam
do
	# echo $file
	echo $(basename $file).peaks.txt
	bedtools coverage -a ../../input/atac_pbmc_5k_nextgem_peaks.bed -b $file > ../../tmp/count_reads_peaks_output/$(basename $file).peaks.txt
done
