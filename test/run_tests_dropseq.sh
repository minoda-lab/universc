#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

# set up cellranger reference
if [[ ! -f test/cellranger_reference/cellranger-tiny-ref/3.0.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA test/cellranger_reference/cellranger-tiny-ref/3.0.0/star/SA
fi
if [[ ! -f test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA
fi

## test drop-seq data
# unzip input data
if [[ -f test/shared/dropseq-test/*fastq.gz ]]; then
    gunzip -f test/shared/dropseq-test/*fastq.gz
fi
# test manual setup
bash launch_universc.sh -t "nadia" --setup
# remove processed files
if [[ -f SRR1873277_S1_L001_R[12]_001.fastq ]]; then
    rm SRR1873277_S1_L001_R[12]_001.fastq
fi
# call on dropseq with files
bash launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/dropseq-test/SRR1873277_Sample1_R1" \
 --read2 "test/shared/dropseq-test/SRR1873277_Sample1_R2" 

#reset test data (files names)
if [[ -f test/shared/dropseq-test/SRR1873277_S1_L001_R[12]_001.fastq ]]; then
    gzip test/shared/dropseq-test/SRR1873277_S1_L001_R[12]_001.fastq
fi
if [[ ! -f test/shared/dropseq-test/SRR1873277_Sample1_R1.fastq.gz ]]; then
    mv test/shared/dropseq-test/SRR1873277_S1_L001_R1_001.fastq.gz test/shared/dropseq-test/SRR1873277_Sample1_R1.fastq.gz
else
    rm test/shared/dropseq-test/SRR1873277_S1_L001_R1_001.fastq.gz
fi
if [[ ! -f test/shared/dropseq-test/SRR1873277_Sample1_R2.fastq.gz ]]; then
    mv test/shared/dropseq-test/SRR1873277_S1_L001_R2_001.fastq.gz test/shared/dropseq-test/SRR1873277_Sample1_R1.fastq.gz
else
    rm test/shared/dropseq-test/SRR1873277_S1_L001_R2_001.fastq.gz
fi

# compress all input files
if [[ -f test/shared/dropseq-test/*fastq ]]; then
    gzip test/shared/dropseq-test/*fastq
fi
