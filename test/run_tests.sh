#!/bin/bash

# run tests in home directory (writeable)
cd ~

# reset barcodes for test
bash /universc/launch_universc.sh -t "10x" --setup

## test 10x data
# unzip input data
if [[ ! -f /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
# test cellranger call
cellranger testrun --id="tiny-test"
# unzip input data
gunzip -fk /universc/test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz
gunzip -fk /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/barcodes/3M-february-2018.txt.gz 
# test cellranger call
cellranger count --id="tiny-count-v3" \
 --fastqs="/cellranger-3.0.2.9001/cellranger-tiny-fastq/3.0.0/" --sample="tinygex" \
 --transcriptome="/cellranger-3.0.2.9001/cellranger-tiny-ref/3.0.0"

cellranger count --id="tiny-count-v2" \
 --fastqs="/cellranger-3.0.2.9001/cellranger-tiny-fastq/1.2.0/" --sample="tinygex" \
 --transcriptome="/cellranger-3.0.2.9001/cellranger-tiny-ref/1.2.0"

# call convert on 10x with multiple lanes
bash /universc/launch_universc.sh --id "test-10x-v3" --technology "10x" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "/universc/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "/universc/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002"

bash /universc/launch_universc.sh --id "test-10x-v2" --technology "10x" \
 --reference "/universc/test/cellranger_reference/cellranger-tiny-ref/1.2.0" \
 --file "/universc/test/shared/cellranger-tiny-fastq/1.2.0/tinygex_S1_L001" \
 "/universc/test/shared/cellranger-tiny-fastq/1.2.0/tinygex_S1_L002"

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
gzip -f /universc/test/shared/cellranger-tiny-fastq/3.0.0/*fastq
gzip -f /universc/test/shared/dropseq-test/*fastq
gzip -f /universc/test/shared/mappa-test/*fastq

