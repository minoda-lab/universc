#!/bin/bash

# run tests in home directory (writeable)
cd ~

## test drop-seq data
# unzip input data
gunzip /universc/test/shared/dropseq-test/*fastq.gz
# test manual setup
bash /universc/launch_universc.sh -t "nadia" --setup
# call on dropseq with files
bash /universc/launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/dropseq-test/SRR1873277_S1_L001_R1_001" \
 --read2 "/universc/test/shared/dropseq-test/SRR1873277_S1_L001_R2_001" 

# compress all input files
gzip -f /universc/test/shared/dropseq-test/*fastq

