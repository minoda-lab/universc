#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd
git pull --ff-only origin $(git branch --show-current) 

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

rm -rf  test/shared/dropseq-test/* test/shared/cellranger-tiny-fastq/* test/shared/mappa-test/  test/shared/icell8-test/ test/shared/indrop-v3-test/ test/shared/smartseq3-test/

bash launch_universc.sh --id "test-sciseq" --technology "sciseq" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/sciseq-v3-test/SRR7827205_S1_R1.fastq.gz" \
 --read2 "test/shared/sciseq-v3-test/SRR7827205_S1_R2.fastq.gz" \
 --per-cell-data --jobmode "local" --localcores 1 
