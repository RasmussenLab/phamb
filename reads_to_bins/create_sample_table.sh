#!/bin/bash 

rm -f sample_table.txt
while read -r id; do
	    fq1="fastq/${id}*R1*.fastq.gz"
	    fq2="fastq/${id}*R2*.fastq.gz"
		  echo ${id} ${fq1} ${fq2} >>sample_table.txt
done < sample.list
