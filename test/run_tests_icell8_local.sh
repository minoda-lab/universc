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

unpigz -k test/shared/icell8-test/72618_KU812*fastq.gz

if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R2_001.fastq*
fi
if [ -d test-icell8-72618-KU812-2-lanes ];then
    rm -rf test-icell8-72618_KU812-2-lanes
fi
bash launch_universc.sh --id "test-icell8-72618_KU812-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_KU812_L001_R1_001.fastq" "test/shared/icell8-test/72618_KU812_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_KU812_L001_R2_001.fastq" "test/shared/icell8-test/72618_KU812_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "local" --localcores 2 --localmem 4 
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R2_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S1_L002_R1_001.fastq ]; then
    rename "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/72618_KU812_S1_L002_R2_001.fastq ]; then
    rename "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R2_001.fastq*
fi
if [ -d test-icell8-72618_KU812-1-lane ]; then
    rm -rf test-icell8-72618_KU812-1-lane
fi
bash launch_universc.sh --id "test-icell8-72618_KU812-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_KU812_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_KU812_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "local" --localcores 2 --localmem 4 
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq*
fi
if [ -f test/shared/icell8-test/KU812_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R2_001.fastq*
fi
