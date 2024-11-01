#! /usr/bin/bash

set -eo pipefail

input=$1
outdir=$2

## get SRA ID per sample
grep "RNAseq" $input | cut -f1 -d "," > file

## download sra file of the samples in parallel
# cat file | parallel --bar -j 8 prefetch {} --output-directory $outdir --progress

## create a new directory for converted fastq
mkdir -p $outdir"fastqs"

# A loop for convert sra files to paired end fastq files
for i in `cat file`;
do
	if [ -f $outdir"fastqs/"$i"_1.fastq.gz" ]
		then
			echo "$i already converted"
	else
			echo "converting $i.."
			fasterq-dump $outdir$i --split-files --skip-technical --threads 8 --outdir $outdir"fastqs"
			echo "compressing $i.."
			gzip $outdir"fastqs/"$i"_1.fastq"
			gzip $outdir"fastqs/"$i"_2.fastq"
	fi
done

rm file
