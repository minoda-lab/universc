#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

if [[ ! -f ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/3.0.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/3.0.0/star/SA
fi
if [[ ! -f ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/1.2.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/1.2.0/star/SA
fi

## test on internal icell8 data
if [ -d test-icell8-default ];then
    rm -rf test-icell8-default
fi
if [ -f test/internal/icell8-test/iCELL8_01_S1_L001_R1_001.fastq.gz ]; then
     gunzip -k test/internal/icell8-test/iCELL8_01_S1_L001_R1_001.fastq.gz
fi
if [ -f test/internal/icell8-test/iCELL8_01_S1_L001_R2_001.fastq.gz ]; then
     gunzip -k test/internal/icell8-test/iCELL8_01_S1_L001_R2_001.fastq.gz
fi
if [ -f test/internal/icell8-test/iCELL8_01_S1_L002_R1_001.fastq.gz ]; then
    gunzip -k test/internal/icell8-test/iCELL8_01_S1_L002_R1_001.fastq.gz
fi
if [ -f test/internal/icell8-test/iCELL8_01_S1_L002_R2_001.fastq.gz ]; then
    gunzip -k test/internal/icell8-test/iCELL8_01_S1_L002_R2_001.fastq.gz
fi
# call on icell8 files with default whitelist
bash launch_universc.sh --id "test-icell8-default" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/internal/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/internal/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/internal/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/internal/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --jobmode "sge"

if [ -d test-icell8-custom ];then
    rm -rf test-icell8-custom
fi
# call on icell8 files with custom whitelist and non-standard file names
bash launch_universc.sh --setup -t "icell8"  --barcodefile "test/internal/icell8-test/BarcodeList.txt"
bash launch_universc.sh --id "test-icell8-custom" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/internal/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/internal/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/internal/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/internal/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --barcodefile "test/internal/icell8-test/BarcodeList.txt" \
 --jobmode "sge"

bash launch_universc.sh --id "test-icell8-test-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/internal/icell8-test/test_L001_R1_001.fastq" "test/internal/icell8-test/test_L002_R1_001.fastq" \
 --read2 "test/internal/icell8-test/test_L001_R2_001.fastq" "test/internal/icell8-test/test_L002_R2_001.fastq" \
 --barcodefile "test/internal/icell8-test/WellList.txt" \
 --jobmode "sge"
bash launch_universc.sh --id "test-icell8-test-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/internal/icell8-test/test_S1_L001_R1_001.fastq" \
 --read2 "test/internal/icell8-test/test_S2_L001_R2_001.fastq" \
 --barcodefile "test/internal/icell8-test/WellList.txt" \
 --jobmode "sge"
