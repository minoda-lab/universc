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

# compress all input files
gzip -f /universc/test/shared/cellranger-tiny-fastq/3.0.0/*fastq
