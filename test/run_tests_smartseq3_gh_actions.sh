#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd
git pull origin master

# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0 test/cellranger_reference/cellranger-tiny-ref/3.0.0
cd test/cellranger_reference/cellranger-tiny-ref/
cellranger mkref --genome=3.0.0 --fasta=genome-3.0.0.fa --genes=genes-3.0.0.gtf
cd ../../..
make -C test/cellranger_reference/cellranger-tiny-ref reference
rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0

rm -rf  test/shared/dropseq-test/* test/shared/cellranger-tiny-fastq/* test/shared/mappa-test/  test/shared/icell8-test/ test/shared/indrop-v3-test/ test/shared/sciseq-v3-test

# compress all input files
if [[ -f test/shared/smartseq3-test/*fastq ]]; then
    gzip test/shared/smartseq3-test/*fastq
fi

# test manual setup
bash launch_universc.sh -t "smartseq" --setup

if [[ -d input4cellranger_test-smartseq3 ]]; then
    rm -rf input4cellranger_test-smartseq3
fi
if [[ -d test-smartseq3 ]]; then
    rm -rf test-smartseq3
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3" --technology "smartseq" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_diySpike_R1" \
 --read2 "test/shared/smartseq3-test/Smartseq3_diySpike_R2" \
 --jobmode "local" --localcores 1 
