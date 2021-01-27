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
    ln $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA test/cellranger_reference/cellranger-tiny-ref/3.0.0/star/SA
fi
if [[ ! -f test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ]]; then
    ln $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA
fi

rm -rf test/shared/dropseq-test/* test/shared/cellranger-tiny-fastq/* test/shared/mappa-test/
rm -rf  test/cellranger_reference/cellranger-tiny-ref/1.2.0
unpigz -f test/shared/indrop-v3-test/*fastq.gz

if [ -d test-indrop-v3 ];then
    rm -rf test-indrop-v3
fi
bash launch_universc.sh --id "test-indrop-v3" --technology "indrops-v3" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/indrop-v3-test/Undetermined_S0_L001_R1_001.fastq" "test/shared/indrop-v3-test/Undetermined_S0_L002_R1_001.fastq" \
 --read2 "test/shared/indrop-v3-test/Undetermined_S0_L001_R2_001.fastq" "test/shared/indrop-v3-test/Undetermined_S0_L002_R2_001.fastq" \
 --jobmode "local" --localcores 1 
