#!/bin/bash

# run tests in home directory (writeable)
#cd ~

#export PATH=/home/tom/local/bin/cellranger-3.1.0:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`


## test icell8 data
# unzip input data
#gzip -fk test/shared/icell8-test/*fastq
gunzip -fk test/shared/icell8-test/*fastq.gz
# call on icell8 files with default whitelist
#bash launch_universc.sh --id "test-icell8-default" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --jobmode "sge"
#restore file names for test job
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/icell8-test/test_FL_R1.fastq
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/icell8-test/test_FL_R2.fastq

# call on icell8 files with custom whitelist and non-standard file names
bash launch_universc.sh --setup -t "icell8"  --barcodefile "test/shared/icell8-test/BarcodeList.txt"
bash launch_universc.sh --id "test-icell8-custom" --technology "iCell8" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/icell8-test/iCELL8_01_S1_L001_R1_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R1_001.fastq" \
 --read2 "test/shared/icell8-test/iCELL8_01_S1_L001_R2_001.fastq" "test/shared/icell8-test/iCELL8_01_S1_L002_R2_001.fastq" \
 --barcodefile "test/shared/icell8-test/BarcodeList.txt" \
 --jobmode "sge"
#restore file names for test job
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/icell8-test/test_FL_R1.fastq
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/icell8-test/test_FL_R2.fastq

# compress all input files
gzip -f test/shared/icell8-test/*fastq

