#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0 test/cellranger_reference/cellranger-tiny-ref/3.0.0
cd test/cellranger_reference/cellranger-tiny-ref/
cellranger mkref --genome=3.0.0 --fasta=genome-3.0.0.fa --genes=genes-3.0.0.gtf
cd ../../..
make -C test/cellranger_reference/cellranger-tiny-ref reference 
rm -rf test/cellranger_reference/cellranger-tiny-ref/3.0.0

rm -rf  test/shared/dropseq-test/* test/shared/mappa-test/  test/shared/icell8-test/ test/shared/smartseq3-test/ test/shared/indrop-v3-test/ test/shared/sciseq-v3-test

# reset barcodes for test
bash launch_universc.sh -t "10x" --setup

if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi

cellranger count --id="tiny-count-v2" \
 --fastqs="test/shared/cellranger-tiny-fastq/1.2.0/" --sample="" --chemistry="threeprime" \
 --transcriptome="test/cellranger_reference/cellranger-tiny-ref/1.2.0" \
 --jobmode "local" --localcores 1 
