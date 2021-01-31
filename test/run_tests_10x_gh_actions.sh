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
rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0

rm -rf  test/shared/dropseq-test/* test/shared/mappa-test/  test/shared/icell8-test/ test/shared/smartseq3-test/ test/shared/indrop-v3-test/

# reset barcodes for test
bash launch_universc.sh -t "10x" --setup

## test 10x data
# unzip input data
if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
if [[ -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt ]]; then
    rm -rf ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi

# call convert on 10x with multiple lanes
if [[ -d test-10x-v3 ]]; then
    rm -rf test-10x-v3
fi
if [[ -d input4cellranger_test-10x-v3 ]]; then
    rm -rf input4cellranger_test-10x-v3
fi
# unzip input data
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001/*fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001/*fastq.gz
fi
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002/*fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002/*fastq.gz
fi
bash launch_universc.sh --id "test-10x-v3" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
  --jobmode "local" --localcores 1 

# compress all input files
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq ]]; then
    gzip -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq
fi
