
#!/bin/bash
bash ~/.bashrc

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd
##git pull --ff-only origin $(git branch --show-current) 

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

while [[ ! -f /home/tom/repos/universc/test-10x-hek293t/outs/metrics_summary.csv  ]]; do
    echo "wait ..."
    sleep 60
done

# test manual setup
bash launch_universc.sh -t "smartseq3" --setup --barcodefile "whitelists/test_bcs_smartseq3_hek293t.txt"

if [[ -d input4cellranger_test-smartseq3-hek293t-test-chr21 ]]; then
    rm -rf input4cellranger_test-smartseq3-hek293t-test-chr21
fi
if [[ -d test-smartseq3-hek293t-test-chr21 ]]; then
    rm -rf test-smartseq3-hek293t-test-chr21
fi

if [[ ! -f whitelists/SmartSeq3_test_barcodes.txt ]]; then
    gunzip -k whitelists/SmartSeq3_test_barcodes.txt.gz
fi


#gzip -fk test/shared/smartseq3-test/Smartseq3_diySpike_[IR][12].fastq

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-hek293t-test-chr21" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --barcodefile "whitelists/test_bcs_smartseq3_hek293t.txt" \
 --read1 "test/shared/smartseq3-test/test-smartseq-hek293t_S2_L001_R1_001.fastq" \
 --read2 "test/shared/smartseq3-test/test-smartseq-hek293t_S2_L001_R2_001.fastq" \
 --index1 "test/shared/smartseq3-test/test-smartseq-hek293t_S2_L001_I1_001.fastq" \
 --index2 "test/shared/smartseq3-test/test-smartseq-hek293t_S2_L001_I2_001.fastq" \
 --per-cell-data --jobmode "local" --verbose #--as-is # "local" --localcores 1
