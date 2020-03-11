# run tests in home directory (writeable)
cd ~

## check cellranger installed
cellranger count ---version

gunzip /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/barcodes/3M-february-2018.txt.gz 


## check convert is installed
bash /universc/launch_universc.sh -v

## test 10x data
# unzip input data
gunzip /universc/test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz
# test cellranger call
cellranger count --id="tiny-test" \
 --fastqs="/cellranger-3.0.2.9001/cellranger-tiny-fastq/3.0.0/" --sample="tinygex" \
 --transcriptome="/cellranger-3.0.2.9001/cellranger-tiny-ref/3.0.0/"

# call convert on 10x with multiple lanes
bash /universc/launch_universc.sh --id "test-10x" --technology "10x" \
 --reference "/universc/test/shared/cellranger-tiny-ref/3.0.0" \
 --file "/universc/test/shared/cellranger-tiny-fastq/3.0.0/tinygex"

## test drop-seq data
# unzip input data
gunzip universc/test/shared/dropseq-test/*fastq.gz
# test manual setup
bash /universc/launch_universc.sh -t "nadia" --setup
# call on dropseq with files
bash /universc/launch_universc.sh --id "test-dropseq" --technology "nadia" \
 --reference "/universc/test/shared/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/dropseq-test/SRR1873277_S1_L001_R1_001" \
 --read2 "/universc/test/shared/dropseq-test/SRR1873277_S1_L001_R2_001" 

## test icell8 data
# unzip input data
gunzip universc/test/shared/mappa-test/*fastq.gz
# call on icell8 files with custom whitelist and non-standard file names
bash /universc/launch_universc.sh --id "test-icell8" --technology "iCell8" \
 --reference "/universc/test/shared/cellranger-tiny-ref/3.0.0" \
 --read1 "/universc/test/shared/mappa-test/test_FL_R1.fastq" \
 --read2 "/universc/test/shared/mappa-test/test_FL_R2.fastq" \
 --barcodefile "/universc/test/shared/dropseq-test/99999_mappa_test_selected_WellList.TXT"

#restore file names for test job
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R1_001.fastq /universc/test/shared/mappa-test/test_FL_R1.fastq
mv /universc/test/shared/mappa-test/test_FL_S1_L001_R2_001.fastq /universc/test/shared/mappa-test/test_FL_R2.fastq

# compress all input files
gzip /universc/test/shared/cellranger-tiny-fastq/3.0.0/*fastq
gzip universc/test/shared/dropseq-test/*fastq
gzip universc/test/shared/mappa-test/*fastq
