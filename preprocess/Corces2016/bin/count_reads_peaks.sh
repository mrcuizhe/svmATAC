#!/bin/bash

#BSUB -J count_reads[1-2034]
#BSUB -o count_reads_peaks.out
#BSUB -e count_reads_peaks.err
#BSUB -We 5
#BSUB -q vshort

bampath=../tmp/corces2016_bam_sorted_nondup
# dirlist=(`ls $bampath*.bam`)
# echo ${dirlist[$LSB_JOBINDEX-1]}
mkdir -p ../tmp/count_reads_peaks_output
for file in $bampath/*
do
	# echo $file
	echo ../tmp/count_reads_peaks_output/$(basename $file).peaks.txt
	bedtools coverage -a ../input/corces2016_cicero.bed -b $file > ../tmp/count_reads_peaks_output/$(basename $file).peaks.txt
done
