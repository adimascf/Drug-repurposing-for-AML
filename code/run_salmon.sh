#! /usr/bin/bash

set -eo pipefail


indir=$1
outdir=$2
index=$3

ls $indir | grep "_1.fastq.gz" | sed 's/_1.fastq.gz//' > file

mkdir -p ${outdir}/quants

for i in `cat file`;
do
	salmon quant -i $index -l A \
         -1 ${indir}${i}_1.fastq.gz \
         -2 ${indir}${i}_2.fastq.gz \
         -p 8 --validateMappings -o ${outdir}/quants/${i}_quant
done
