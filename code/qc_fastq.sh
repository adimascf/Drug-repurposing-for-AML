#! /usr/bin/bash

set -eo pipefail

indir=$1
outdir=$2
threads=$3


ls $1 | grep "_1.fastq.gz" | sed 's/_1.fastq.gz//' > file

mkdir -p $outdir/fastqs-qcd/
mkdir -p $outdir/json

cat file | parallel --bar -j $threads \
        fastp -i $indir/{}_1.fastq.gz -I $indir/{}_2.fastq.gz \
                -o $outdir/fastqs-qcd/{}_1.fastq.gz -O $outdir/fastqs-qcd/{}_2.fastq.gz \
                --json $outdir/json/{}.fastp.json --html $outdir/json/{}.fastp.html \
	       	--correction --thread 6 --length_required 75

rm file
