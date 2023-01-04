#!/bin/bash

cat sample_sheet* | sed '1d' > NoHead.sampleSheet.csv;

echo >> NoHead.sampleSheet.csv;

SAMPLEFILE=NoHead.sampleSheet.csv;

FILEout=sampleSheet.csv;

mkdir -p "$PWD"/RAWREADS;

NanoRawDir="$PWD"/fastq_pass;

RAW="$PWD"/RAWREADS;

SAMPLElist="$PWD"/samples.txt;

##  count for walking through the input file
count1=1;
countS=$(cat "$SAMPLEFILE" | wc -l);

	while [ "$count1" -le "$countS" ];do
##  select the line for parsing
		LINE=$(cat "$SAMPLEFILE" | awk 'NR=='"$count1");
##  variables needed for downstream analysis data
			sample=$(echo "$LINE" | cut -f9 -d',');
			barcode=$(echo "$LINE" | cut -f8 -d',');
##  concatenate all teh individual raw output files into 1 file with the name of the sample
			cat "$NanoRawDir"/"$barcode"/*.gz > "$RAW"/"$sample".fastq.gz;
##  create a sample file for downstream processing
			echo "$sample" >> "$SAMPLElist";
##  counter
	count1=$((count1+1));
	done

exit 1
