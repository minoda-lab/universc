
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

# compress all input files
if [[ -f test/shared/smartseq3-test/Smartseq3.diySpike.R1.fastq ]]; then
echo    gzip test/shared/smartseq3-test//Smartseq3.diySpike*fastq
fi

rename -f "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_[IR][12]_001.fastq
rename -f "s/_001//" test/shared/smartseq3-test/Smartseq3_diySpike_[IR][12]_001.fastq

#rename "s/Smartseq3.diySpike./Smartseq3_diySpike_/g"  test/shared/smartseq3-test/Smartseq3.diySpike.[IR][12].fastq
#rename "s/Smartseq3.diySpike./Smartseq3_diySpike_/g"  test/shared/smartseq3-test/Smartseq3.diySpike.L001.[IR][12].fastq

while [[ ! -f /home/tom/repos/universc/test-10x-hek293t/outs/metrics_summary.csv  ]]; do
    echo "wait ..."
    sleep 60
done

# test manual setup
bash launch_universc.sh -t "smartseq3" --setup --barcodefile "whitelists/test_bcs_small.txt"

if [[ -d input4cellranger_test-smartseq3-diyspike-hg38 ]]; then
    rm -rf input4cellranger_test-smartseq3-diyspike-hg38
fi
if [[ -d test-smartseq3-diyspike-hg38 ]]; then
    rm -rf test-smartseq3-diyspike-hg38
fi

if [[ ! -f whitelists/SmartSeq3_test_barcodes.txt ]]; then
    gunzip -k whitelists/SmartSeq3_test_barcodes.txt.gz
fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-diyspike-hg38" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq.gz" \
 --read2 "test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R2_001.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_I1_001.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_I2_001.fastq" \
 --barcodefile "whitelists/Smartseq3_diySpike_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_diySpike_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_diySpike_[IR][12]_001.fastq*
fi
