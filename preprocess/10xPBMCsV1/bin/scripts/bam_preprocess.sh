#! /bin/bash

module load SAMtools/1.9-foss-2018b
mkdir -p ../../tmp/10xpbmc_bam_sorted_nondup
mkdir -p ../../tmp/10xpbmc_bam_sorted
for file in ../../tmp/sc_bam_dedup/*.bam
do
	bamName=`echo $file|cut -f5 -d"/" |cut -f1 -d"."`
	if [[ ! -f ../../tmp/10xpbmc_bam_sorted_nondup/"$bamName"_sorted_nondup.bam ]]
	then
		#sort bam
		samtools sort -@ 8 $file > ../../tmp/10xpbmc_bam_sorted/"$bamName"_sorted.bam
		#remove duplicates
		java -jar picard.jar MarkDuplicates I=../../tmp/10xpbmc_bam_sorted/"$bamName"_sorted.bam O=../../tmp/10xpbmc_bam_sorted_nondup/"$bamName"_sorted_nondup.bam M=../../tmp/10xpbmc_bam_sorted_nondup/"$bamName"_dup_metrics.txt REMOVE_DUPLICATES=true TMP_DIR=../../tmp/tmp_markduplicates USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
		echo "Remove duplicate for "$bamName" is finished"
	fi
done
