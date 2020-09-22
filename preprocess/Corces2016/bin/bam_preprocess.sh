#! /bin/bash

mkdir -p ../tmp/corces2016_bam_filtered_sorted_nondup
mkdir -p ../tmp/corces2016_bam_sorted_nondup
for file in ../tmp/corces2016_bam/*.bam
do
	bamName=`echo $file|cut -f4 -d"/" |cut -f1 -d"."`
	if [[ ! -f ../tmp/corces2016_bam_filtered_sorted_nondup/"$bamName"_sorted_nondup_filtered.bam ]]
	then
		#sort bam
		samtools sort -@ 8 ../tmp/corces2016_bam/"$bamName".bam -o ../tmp/corces2016_bam_sorted_nondup/"$bamName"_sorted.bam 
		echo "Sorting for "$bamName" is finished"
		#remove duplicates
		java -jar ../bin/picard.jar MarkDuplicates I=../tmp/corces2016_bam_sorted_nondup/"$bamName"_sorted.bam O=../tmp/corces2016_bam_sorted_nondup/"$bamName"_sorted_nondup.bam M=../tmp/corces2016_bam_sorted_nondup/"$bamName"_dup_metrics.txt REMOVE_DUPLICATES=true TMP_DIR=../tmp/tmp_markduplicates USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
		rm ../tmp/corces2016_bam_sorted_nondup/"$bamName"_sorted.bam
		echo "Remove duplicate for "$bamName" is finished"
		#two ends mapped in proper pair && MAPQ>=30 && fragment length<=1000
		samtools view -h -S -F 0x0002 -q 30 ../tmp/corces2016_bam_sorted_nondup/"$bamName"_sorted_nondup.bam -b > ../tmp/corces2016_bam_filtered_sorted_nondup/"$bamName"_sorted_nondup_filtered.bam
		echo "Filtering for "$bamName" is finished"
	fi
done