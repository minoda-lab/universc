#!/bin/bash

# run tests in home directory (writeable)
#cd ~

#export PATH=/home/tom/local/bin/cellranger-3.1.0:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

## test drop-seq data
# unzip input data
gunzip -fk test/shared/dropseq-test/*fastq.gz
# test manual setup
bash launch_universc.sh -t "nadia" --setup
# call on dropseq with files
bash launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/dropseq-test/SRR1873277_S1_L001_R1_001" \
 --read2 "test/shared/dropseq-test/SRR1873277_S1_L001_R2_001" 

# compress all input files
gzip -f test/shared/dropseq-test/*fastq

