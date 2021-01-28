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

rm -rf test/shared/cellranger-tiny-fastq/* test/shared/mappa-test/ test/shared/icell8-test/ test/shared/indrop-v3-test/
rm -rf  test/cellranger_reference/cellranger-tiny-ref/1.2.0

# compress all input files
if [[ -f test/shared/smartseq3-test/*fastq ]]; then
    gzip test/shared/smartseq3-test/*fastq
fi

# reset test data (compress)
 if [[ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq ]]; then
   gzip -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq
fi
if [[ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq ]]; then
    gzip -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq
fi
# reset test data (files names)
if [[ ! -f test/shared/smartseq3-test/Smartseq3_diySpike_R1.fastq.gz ]]; then
    mv test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz test/shared/smartseq3-test/Smartseq3_diySpike_R1.fastq.gz
else
    rm test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz
fi
if [[ ! -f test/shared/smartseq3-test/Smartseq3_diySpike_R2.fastq.gz ]]; then
    mv test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq.gz test/shared/smartseq3-test/Smartseq3_diySpike_R2.fastq.gz
else
    rm test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq.gz
fi

## test drop-seq data
# unzip input data
if [[ -f test/shared/smartseq3-test/*fastq.gz ]]; then
    unpigz -f test/shared/smartseq3-test/*fastq.gz
fi
# test manual setup
bash launch_universc.sh -t "smartseq" --setup

# remove processed files
if [[ -f SRR1873277_S1_L001_R[12]_001.fastq ]]; then
    rm SRR1873277_S1_L001_R[12]_001.fastq
fi
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

# reset test data (compress)
 if [[ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq ]]; then
   gzip -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq
fi
if [[ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq ]]; then
    gzip -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq
fi
# reset test data (files names)
if [[ ! -f test/shared/smartseq3-test/Smartseq3_diySpike_R1.fastq.gz ]]; then
    mv test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz test/shared/smartseq3-test/Smartseq3_diySpike_R1.fastq.gz
else
    rm test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz
fi
if [[ ! -f test/shared/smartseq3-test/Smartseq3_diySpike_R2.fastq.gz ]]; then
    mv test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq.gz test/shared/smartseq3-test/Smartseq3_diySpike_R2.fastq.gz
else
    rm test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq.gz
fi

# compress all input files
if [[ -f test/shared/smartseq3-test/*fastq ]]; then
    gzip test/shared/smartseq3-test/*fastq
fi
