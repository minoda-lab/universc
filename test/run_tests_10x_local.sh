#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

# used to export to PATH for testing on SGE server
export PATH=${HOME}/local/bin/cellranger-2.1.0:$PATH

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`

# set up cellranger reference
if [[ ! -f test/cellranger_reference/cellranger-tiny-ref/3.0.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA test/cellranger_reference/cellranger-tiny-ref/3.0.0/star/SA
fi
if [[ ! -f test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ]]; then
    rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA test/cellranger_reference/cellranger-tiny-ref/1.2.0/star/SA
fi

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

# reset barcodes for test
bash launch_universc.sh -t "10x" --setup

## test 10x data
# unzip input data
ls ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt*
if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
rm -rf ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
if [[ -d tiny-test ]]; then
    rm -rf tiny-test
fi
# test cellranger call
if [[ $cellrangerversion == "3.0.2" ]]; then
    cellranger testrun --id="tiny-test"
fi
# unzip input data
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz ]]; then
    gunzip -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz
fi
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq ]]; then
    rm -rf test/shared/cellranger-tiny-fastq/3.0.0/*fastq.gz
fi
# test cellranger call
if [[ -d tiny-count-v3 ]]; then
    rm -rf tiny-count-v3
fi
if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
if [[ -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt ]]; then
    rm -rf ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
cellranger count --id="tiny-count-v3" \
 --fastqs="test/shared/cellranger-tiny-fastq/3.0.0/" --sample="tinygex" \
 --transcriptome="test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --jobmode "local" --localcores 2 --localmem 4 
if [[ -d tiny-count-v2 ]]; then
    rm -rf tiny-count-v2
fi 
if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz ]]; then
    gzip -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
if [[ -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt ]]; then
    rm -rf ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/3M-february-2018.txt
fi
cellranger count --id="tiny-count-v2" \
 --fastqs="test/shared/cellranger-tiny-fastq/1.2.0/" --sample="" --chemistry="threeprime" \
 --transcriptome="test/cellranger_reference/cellranger-tiny-ref/1.2.0" \
 --jobmode "local" --localcores 2 --localmem 4 

# call convert on 10x with multiple lanes
if [[ -d test-10x-v3 ]]; then
    rm -rf test-10x-v3
fi
if [[ -d input4cellranger_test-10x-v3 ]]; then
    rm -rf input4cellranger_test-10x-v3
fi
# unzip input data
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq.gz
    rm -rf test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq.gz
fi
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq.gz
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq.gz
fi
bash launch_universc.sh --id "test-10x-v3" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
 --per-cell-data --jobmode "local" --localcores 2 --localmem 4

# call convert on 10x with multiple lanes with compressed files
if [[ -d test-10x-v3-zip ]]; then
    rm -rf test-10x-v3-zip
fi
if [[ -d input4cellranger_test-10x-v3-zip ]]; then
    rm -rf input4cellranger_test-10x-v3-zip
fi
# unzip input data
if [[ ! -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz ]]; then
    pigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq
fi
if [[ ! -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz ]]; then
    pigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq
    rm -rf test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq
fi
bash launch_universc.sh --id "test-10x-v3-zip" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
 --per-cell-data --jobmode "local" --localcores 2 --localmem 4

# call convert on 10x with multiple lanes with *some* compressed files
if [[ -d test-10x-v3-zip-mix ]]; then
    rm -rf test-10x-v3-zip-mix
fi
if [[ -d input4cellranger_test-10x-v3-zip-mix ]]; then
    rm -rf input4cellranger_test-10x-v3-zip-mix
fi
# unzip input data
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq.gz
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001*fastq.gz
fi
if [[ ! -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz ]]; then
    pigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002*fastq
fi
bash launch_universc.sh --id "test-10x-v3-zip-mix" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
 --per-cell-data --jobmode "local" --localcores 2 --localmem 4

# call convert on 10x with multiple lanes with *some* compressed files
if [[ -d test-10x-v3-zip-mix-2 ]]; then
    rm -rf test-10x-v3-zip-mix-2
fi
if [[ -d input4cellranger_test-10x-v3-zip-mix-2 ]]; then
    rm -rf input4cellranger_test-10x-v3-zip-mix-2
fi
# unzip input data
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq.gz
    pigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R2_001.fastq
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R2_001.fastq
fi
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz ]]; then
    unpigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq.gz
    pigz -f test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R2_001.fastq
    rm -rf  test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R2_001.fastq
fi
bash launch_universc.sh --id "test-10x-v3-zip-mix-2" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
 --per-cell-data --jobmode "local" --localcores 2 --localmem 4

# call convert on 10x with multiple lanes with *some* compressed files
if [[ -d test-10x-v3-zip-mix-3 ]]; then
    rm -rf test-10x-v3-zip-mix-3
fi
if [[ -d input4cellranger_test-10x-v3-zip-mix-3 ]]; then
    rm -rf input4cellranger_test-10x-v3-zip-mix-3
fi
bash launch_universc.sh --id "test-10x-v3-zip-mix-3" --technology "10x" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001" \
 "test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002" \
 --per-cell-data --jobmode "local" --localcores 2 --localmem 4

# compress all input files
if [[ -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq ]]; then
    gzip -f test/shared/cellranger-tiny-fastq/3.0.0/*fastq
fi
