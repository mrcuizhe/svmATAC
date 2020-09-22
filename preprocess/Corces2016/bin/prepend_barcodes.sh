#!/bin/bash
OUTPUT='../tmp/sc-bams_barcodes'
mkdir -p $OUTPUT
FILEPATH='../tmp/corces2016_bam_sorted_nondup/'
METAFILE='../tmp/corces2016_barcode_metadata.tsv'
{
	read
	while IFS=$'\t' read -r -a myArray
	do
		echo "${myArray[0]}"
		echo "${myArray[1]}"
		# echo $FILEPATH${myArray[1]}.st.bam
		# echo $OUTPUT/${myArray[1]}.with_barcode.bam
	./prependBarcode $FILEPATH${myArray[1]}_sorted_nondup.bam ${myArray[0]} $OUTPUT/${myArray[1]}.with_barcode.bam
	done
}< $METAFILE



