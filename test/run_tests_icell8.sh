#!/bin/bash

# run tests in home directory (writeable)
#cd ~
# run tests in universc directory
#cd $(dirname ${BASH_SOURCE[0]})/..
pwd
# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-3.0.2:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

## test icell8 data
# unzip input data
#gzip -k test/shared/icell8-test/*fastq
#gzip -k test/shared/icell8-test/*/*/*fastq 
#gunzip -k test/shared/icell8-test/*fastq
#gunzip -k test/shared/icell8-test/*/*/*fastq 
#if [[ -f test/shared/icell8-test/*fastq.gz ]]; then
#    gunzip -fk test/shared/icell8-test/*fastq.gz
#fi
if [ -d test-icell8-default ];then
    rm -rf test-icell8-default
fi
# call on icell8 files with default whitelist
bash launch_universc.sh --id "test-icell8-default" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --jobmode "sge"

if [ -d test-icell8-custom ];then
    rm -rf test-icell8-custom
fi
# call on icell8 files with custom whitelist and non-standard file names
bash launch_universc.sh --setup -t "icell8"  --barcodefile "test/shared/icell8-test/BarcodeList.txt"
bash launch_universc.sh --id "test-icell8-custom" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/BarcodeList.txt" \
 --jobmode "sge"

if [ -f test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_A375_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_A375_S2_L002_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_A375_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_A375_S2_L002_R2_001.fastq
fi
if [ -d test-icell8-72618-A375-2-lanes ];then
    rm -rf test-icell8-72618-A375-2-lanes
fi
bash launch_universc.sh --id "test-icell8-72618_A375-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_A375_L001_R1_001.fastq" "test/shared/icell8-test/72618_A375_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_A375_L001_R2_001.fastq" "test/shared/icell8-test/72618_A375_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
if [ -f test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_A375_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/A375_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_A375_S1_L001_R2_001.fastq
fi
if [ -d test-icell8-A375-1-lane ];then
    rm -rf test-icell8-A375-1-lane
fi
bash launch_universc.sh --id "test-icell8-72618_A375-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_A375_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_A375_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"

if [ -f test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_HCT116_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_HCT116_S2_L002_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_HCT116_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_HCT116_S2_L002_R2_001.fastq
fi
if [ -d test-icell8-72618-HCT116-2-lanes ];then
    rm -rf test-icell8-72618-HCT116-2-lanes
fi
bash launch_universc.sh --id "test-icell8-72618_HCT116-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_HCT116_L001_R1_001.fastq" "test/shared/icell8-test/72618_HCT116_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_HCT116_L001_R2_001.fastq" "test/shared/icell8-test/72618_HCT116_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
if [ -f test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_HCT116_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/HCT116_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_HCT116_S1_L001_R2_001.fastq
fi
if [ -d test-icell8-HCT116-1-lane ];then
    rm -rf test-icell8-HCT116-1-lane
fi
bash launch_universc.sh --id "test-icell8-72618_HCT116-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_HCT116_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_HCT116_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"

if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_KU812_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_KU812_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_KU812_S2_L002_R2_001.fastq
fi
if [ -d test-icell8-72618-KU812-2-lanes ];then
    rm -rf test-icell8-72618-KU812-2-lanes
fi
bash launch_universc.sh --id "test-icell8-72618_KU812-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_KU812_L001_R1_001.fastq" "test/shared/icell8-test/72618_KU812_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_KU812_L001_R2_001.fastq" "test/shared/icell8-test/72618_KU812_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
if [ -f test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/KU812_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_KU812_S1_L001_R2_001.fastq
fi
if [ -d test-icell8-KU812-1-lane ];then
    rm -rf test-icell8-KU812-1-lane
fi
bash launch_universc.sh --id "test-icell8-72618_KU812-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_KU812_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_KU812_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"

if [ -f test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_NCI-H2452_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_NCI-H2452_S2_L002_R1_001.fastq
fi
if [ -f test/shared/icell8-test/72618_NCI-H2452_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/72618_NCI-H2452_S2_L002_R2_001.fastq
fi
if [ -d test-icell8-72618-NCI-H2452-2-lanes ];then
    rm -rf test-icell8-72618-NCI-H2452-2-lanes
fi
bash launch_universc.sh --id "test-icell8-72618_NCI-H2452-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_NCI-H2452_L001_R1_001.fastq" "test/shared/icell8-test/72618_NCI-H2452_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_NCI-H2452_L001_R2_001.fastq" "test/shared/icell8-test/72618_NCI-H2452_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
if [ -f test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_NCI-H2452_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/NCI-H2452_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/72618_NCI-H2452_S1_L001_R2_001.fastq
fi
if [ -d test-icell8-NCI-H2452-1-lane ];then
    rm -rf test-icell8-NCI-H2452-1-lane
fi
bash launch_universc.sh --id "test-icell8-72618_NCI-H2452-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/72618_NCI-H2452_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/72618_NCI-H2452_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"

if [ -f test/shared/icell8-test/test_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/test_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/test_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/test_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/test_S2_L002_R1_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/test_S2_L002_R1_001.fastq
fi
if [ -f test/shared/icell8-test/test_S2_L002_R2_001.fastq ]; then
    rename -n "s/_S2_L002/_L002/" test/shared/icell8-test/test_S2_L002_R2_001.fastq
fi
if [ -d test-icell8-72618-$-2-lanes ];then
    rm -rf test-icell8-72618-$-2-lanes
fi
bash launch_universc.sh --id "test-icell8-test-2-lanes" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/test_L001_R1_001.fastq" "test/shared/icell8-test/test_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/test_L001_R2_001.fastq" "test/shared/icell8-test/test_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
if [ -f test/shared/icell8-test/test_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/test_S1_L001_R1_001.fastq
fi
if [ -f test/shared/icell8-test/$_S1_L001_R2_001.fastq ]; then
    rename "s/_S1_L001/_L001/" test/shared/icell8-test/test_S1_L001_R2_001.fastq
fi
if [ -d test-icell8-$-1-lane ];then
    rm -rf test-icell8-$-1-lane
fi
bash launch_universc.sh --id "test-icell8-test-1-lane" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/test_L001_R1_001.fastq" \
 --read2 "test/shared/icell8-test/test_L001_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/WellList.txt" \
 --jobmode "sge"
