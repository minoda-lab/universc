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

# compress all input files
if [[ -f test/shared/smartseq3-test/*fastq ]]; then
    gzip test/shared/smartseq3-test/*fastq
fi

mv  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_L001_R1_001.fastq.gz  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_R1.fastq.gz
mv  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_L001_R2_001.fastq.gz  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_R2.fastq.gz
mv  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_L001_I1_001.fastq.gz  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_I1.fastq.gz
mv  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_L001_I2_001.fastq.gz  test/shared/smartseq3-test/Smartseq3_diySpike5_S1_I2.fastq.gz

# test manual setup
bash launch_universc.sh -t "smartseq3" --setup

if [[ -d input4cellranger_test-smartseq3 ]]; then
    rm -rf input4cellranger_test-smartseq3
fi
if [[ -d test-smartseq3 ]]; then
    rm -rf test-smartseq3
fi

if [[ ! -f ${whitelistdir}/SmartSeq3_test_barcodes.txt ]]; then
    gunzip -k ${whitelistdir}/SmartSeq3_test_barcodes.txt.gz
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3" --technology "smartseq3" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_diySpike_R1" \
 --read2 "test/shared/smartseq3-test/Smartseq3_diySpike_R2" \
 --barcodefile "${whitelistdir}/SmartSeq3_test_barcodes.txt" \
 --per-cell-data --jobmode "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz ]; then
    rename "s/_S1_L001/_S1/" test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_diySpike_S1_[IR][12]_001.fastq*
fi
