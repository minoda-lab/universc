#!/bin/bash

# run tests in home directory (writeable)
cd ~

## test icell8 data
# unzip input data
gunzip -fk /universc/test/shared/mappa-test/*fastq.gz
# call on icell8 files with default whitelist
bash /universc/launch_universc.sh --id "test-icell8" --technology "iCell8" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/mappa-test/test_FL_R1.fastq" \
 --read2 "/universc/test/shared/mappa-test/test_FL_R2.fastq" \
 --chemistry "SC3Pv2"
#restore file names for test job
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/mappa-test/test_FL_R1.fastq
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/mappa-test/test_FL_R2.fastq

# call on icell8 files with custom whitelist and non-standard file names
bash /universc/launch_universc.sh --id "test-icell8-custom" --technology "iCell8" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/mappa-test/test_FL_R1.fastq" \
 --read2 "/universc/test/shared/mappa-test/test_FL_R2.fastq" \
 --barcodefile "/universc/test/shared/mappa-test/mappa_barcodes.txt"  \
 --chemistry "SC3Pv2"
#restore file names for test job
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/mappa-test/test_FL_R1.fastq
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/mappa-test/test_FL_R2.fastq

# compress all input files
gzip -f /universc/test/shared/mappa-test/*fastq

