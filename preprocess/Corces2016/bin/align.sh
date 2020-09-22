#!/bin/bash
# Index Reference Genome
snaptools index-genome	\
	--input-fasta=../input/hg19.fa	\
	--output-prefix=hg19	\
   --aligner=bwa	\
	--num-threads=8

# Align
cd ../tmp/corces2016_bam
for sra in ` less ../input/sraList2.txt `
do
	echo $sra
	snaptools align-paired-end	\
		--input-reference=/Volumes/LACIE_SHARE/project/2020-05-27-svmPipeline_generate_corces2016-snap-full_AND_run_cicero_co_accessbility/input/hg19.fa	\
		--input-fastq1=../../input/fq_data/"$sra"_1.fastq	\
		--input-fastq2=../../input/fq_data/"$sra"_2.fastq	\
		--output-bam="$sra".bam	\
		--aligner=bwa	\
		--read-fastq-command=cat	\
		--min-cov=0	\
		--num-threads=8	\
		--if-sort=True	\
		--tmp-folder=../	\
		--overwrite=TRUE    
done
