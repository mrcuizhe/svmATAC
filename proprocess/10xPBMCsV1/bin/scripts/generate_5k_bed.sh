#!/bin/bash
while read line
do
	chrom=`echo $line|cut -f1 -d" "`
	maxLen=`echo $line|cut -f2 -d" "`
	length=1
	while [[ "$length" -le "$maxLen" ]]
	do
		let "length2=length+4999"
    	echo $chrom"\t"$length"\t"$length2 >> ../../input/corces2016_5k_full.bed
    	let "length=length+5000"
    done
done < ../../input/hg19/hg19.chrom.sizes