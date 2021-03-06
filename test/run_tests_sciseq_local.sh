#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd
##git pull --ff-only origin $(git branch --show-current) 

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

##unpigz -k test/shared/sciseq-v3-test/SRR7827205*fastq.gz

mv  test/shared/sciseq-v3-test/SRR7827205_S1_L001_R1_001.fastq.gz  test/shared/sciseq-v3-test/SRR7827205_S1_R1.fastq.gz
mv  test/shared/sciseq-v3-test/SRR7827205_S1_L001_R2_001.fastq.gz  test/shared/sciseq-v3-test/SRR7827205_S1_R2.fastq.gz
mv  test/shared/sciseq-v3-test/SRR7827205_S1_L001_I1_001.fastq.gz  test/shared/sciseq-v3-test/SRR7827205_S1_I1.fastq.gz
mv  test/shared/sciseq-v3-test/SRR7827205_S1_L001_I2_001.fastq.gz  test/shared/sciseq-v3-test/SRR7827205_S1_I2.fastq.gz

#if [ -f test/shared/sciseq-v3-test/SRR7827205_S1_L001_R1_001.fastq ]; then
#    rename "s/_S1_L001/_S1/" test/shared/sciseq-v3-test/SRR7827205_S1_L001_R1_001.fastq*
#    rename "s/_001//" test/shared/sciseq-v3-test/SRR7827205_S1_R1_001.fastq*
#fi
#if [ -f test/shared/sciseq-v3-test/SRR7827205_S1_L001_R2_001.fastq ]; then
#    rename "s/_S1_L001/_S1/" test/shared/sciseq-v3-test/SRR7827205_S1_L001_R2_001.fastq*
#    rename "s/_001//" test/shared/sciseq-v3-test/SRR7827205_S1_R2_001.fastq*
#fi
if [ -d test-sciseq ];then
    rm -rf test-sciseq
fi

if [[ ! -f whitelists/sciseq3_barcode_test.txt ]]; then
    bash launch_universc.sh --id "test-sciseq" --technology "sciseq3" --setup --verbose
    grep "CCGAATCCGACTCCATCGA" whitelists/sciseq3_barcode.txt > whitelists/sciseq3_barcode_test.txt
fi

bash launch_universc.sh --id "test-sciseq" --technology "sciseq" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/sciseq-v3-test/SRR7827205_S1_R1.fastq.gz" \
 --read2 "test/shared/sciseq-v3-test/SRR7827205_S1_R2.fastq.gz" \
 --barcodefile whitelists/sciseq3_barcode_test.txt \
 --per-cell-data --jobmode "local" --localcores 1 --verbose

if [ -f test/shared/sciseq-v3-test/SRR7827205_S1_L001_R1_001.fastq.gz ]; then
    rename "s/_S1_L001/_S1/" test/shared/sciseq-v3-test/SRR7827205_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/sciseq-v3-test/SRR7827205_S1_[IR][12]_001.fastq*
fi
gzip test/shared/sciseq-v3-test/SRR7827205_S1_[IR][12].fastq
