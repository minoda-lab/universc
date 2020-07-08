#!/bin/bash

# run tests in home directory (writeable)
#cd ~
# run tests in universc directory
 cd $(dirname ${BASH_SOURCE[0]})/..

# used to export to PATH for testing on SGE server
#export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

## test drop-seq data
# unzip input data
if [[ -f test/shared/dropseq-test/*fastq.gz ]]; then
    gunzip -f test/shared/dropseq-test/*fastq.gz
fi
# test manual setup
bash launch_universc.sh -t "nadia" --setup
# call on dropseq with files
bash launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/dropseq-test/SRR1873277_S1_L001_R1_001" \
 --read2 "test/shared/dropseq-test/SRR1873277_S1_L001_R2_001" 

# compress all input files
if [[ -f test/shared/dropseq-test/*fastq ]]; then
    gzip -f test/shared/dropseq-test/*fastq
fi
