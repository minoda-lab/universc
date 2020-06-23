#!/bin/bash

# run tests in home directory (writeable)
#cd ~

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

# reset barcodes for test
bash launch_universc.sh -t "10x" --setup

## test 10x data
# unzip input data
if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
# test cellranger call
cellranger testrun --id="tiny-test"
# unzip input data
gunzip -fk test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz
gunzip -fk ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz 
# test cellranger call
cellranger count --id="tiny-count-v3" \
 --fastqs="test/shared/cellranger-tiny-fastq/3.0.0/" --sample="tinygex" \
 --transcriptome="test/cellranger_reference/cellranger-tiny-ref/3.0.0"

cellranger count --id="tiny-count-v2" \
 --fastqs="test/shared/cellranger-tiny-fastq/1.2.0/" --sample="tinygex" \
 --transcriptome="test/cellranger_reference/cellranger-tiny-ref/1.2.0"

# call convert on 10x with multiple lanes
bash launch_universc.sh --id "test-10x-v3" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002"

bash launch_universc.sh --id "test-10x-v2" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/1.2.0" \
 --file "test/shared/cellranger-tiny-fastq/1.2.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/1.2.0/tinygex_S1_L002"

# compress all input files
gzip -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq
