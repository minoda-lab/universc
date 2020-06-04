#!/bin/bash

# run tests in home directory (writeable)
cd ~

## test icell8 data
# unzip input data
gunzip -fk /universc/test/shared/icell8-test/*fastq.gz
# call on icell8 files with default whitelist
bash /universc/launch_universc.sh --id "test-icell8-default" --technology "iCell8" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/icell8-test/iCELL8_01_S1_L001_R1_001 /universc/test/shared/icell8-test/iCELL8_01_S1_L002_R1_001" \
 --read2 "/universc/test/shared/icell8-test/iCELL8_01_S1_L001_R2_001 /universc/test/shared/icell8-test/iCELL8_01_S1_L002_R2_001" 
#restore file names for test job
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/icell8-test/test_FL_R1.fastq
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/icell8-test/test_FL_R2.fastq

# call on icell8 files with custom whitelist and non-standard file names
bash /universc/launch_universc.sh --id "test-icell8-custom" --technology "iCell8" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/icell8-test/iCELL8_01_S1_L001_R1_001 /universc/test/shared/icell8-test/iCELL8_01_S1_L002_R1_001" \
 --read2 "/universc/test/shared/icell8-test/iCELL8_01_S1_L001_R2_001 /universc/test/shared/icell8-test/iCELL8_01_S1_L002_R2_001" 
 --barcodefile "/universc/test/shared/icell8-test/BarcodeList.txt"
#restore file names for test job
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/icell8-test/test_FL_R1.fastq
#mv /universc/test/shared/icell8-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/icell8-test/test_FL_R2.fastq

# compress all input files
gzip -f /universc/test/shared/icell8-test/*fastq

