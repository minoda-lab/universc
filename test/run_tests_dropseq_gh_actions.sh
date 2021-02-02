#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

#rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0 test/cellranger_reference/cellranger-tiny-ref/3.0.0
#cd test/cellranger_reference/cellranger-tiny-ref/
#cellranger mkref --genome=3.0.0 --fasta=genome-3.0.0.fa --genes=genes-3.0.0.gtf
#cd ../../..
#make -C test/cellranger_reference/cellranger-tiny-ref reference 
#rm -rf test/cellranger_reference/cellranger-tiny-ref/1.2.0

#rm -rf test/shared/dropseq-test/* test/shared/cellranger-tiny-fastq/* test/shared/mappa-test/ test/shared/smartseq3-test/ test/shared/indrop-v3-test/ test/shared/sciseq-v3-test

# test manual setup
bash launch_universc.sh -t "nadia" --setup

if [[ -d input4cellranger_test-dropseq ]]; then
    rm -rf input4cellranger_test-dropseq
fi
if [[ -d test-dropseq ]]; then
    rm -rf test-dropseq
fi

# call on dropseq with files
bash launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/dropseq-test/SRR1873277_Sample1_R1" \
 --read2 "test/shared/dropseq-test/SRR1873277_Sample1_R2" \
 --jobmode "local" --localcores 1 
