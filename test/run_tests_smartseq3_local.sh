
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
if [[ -f test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.R1.fastq ]]; then
echo    gzip test/shared/smartseq3-test//Smartseq3.Fibroblasts.GelCut*fastq
fi

#mv  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.L001_R1_001.fastq.gz  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.R1.fastq.gz
#mv  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.L001_R2_001.fastq.gz  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.R2.fastq.gz
#mv  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.L001_I1_001.fastq.gz  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.I1.fastq.gz
#mv  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.L001_I2_001.fastq.gz  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.I2.fastq.gz

rename -f "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq
rename -f "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq

#rename "s/Smartseq3.Fibroblasts.GelCut./Smartseq3_Fibroblasts_GelCut_/g"  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.[IR][12].fastq
#rename "s/Smartseq3.Fibroblasts.GelCut./Smartseq3_Fibroblasts_GelCut_/g"  test/shared/smartseq3-test/Smartseq3.Fibroblasts.GelCut.L001.[IR][12].fastq

# test manual setup
bash launch_universc.sh -t "smartseq3" --setup --barcodefile "whitelists/test_bcs_small.txt"

if [[ false ]]; then

if [[ -d input4cellranger_test-smartseq3-gelcut-small-hard ]]; then
    rm -rf input4cellranger_test-smartseq3-gelcut-small-hard
fi
if [[ -d test-smartseq3-gelcut-small-hard ]]; then
    rm -rf test-smartseq3-gelcut-small-hard
fi

if [[ ! -f whitelists/SmartSeq3_test_barcodes.txt ]]; then
    gunzip -k whitelists/SmartSeq3_test_barcodes.txt.gz
fi


skip="false"
#while [[ ! -f /home/tom/repos/universc/test-smartseq3-gelcut/outs/metrics_summary.csv  ]]; do
#    echo "wait ..."
#    sleep 60
#done

#gzip -fk test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_[IR][12].fastq

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard-5pr2-hg38" --technology "smartseq3" \
 --chemistry "SC5P-R2" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard-5pr1-hg38" --technology "smartseq3" \
 --chemistry "SC5P-R1" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard-hg38" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard-5pr2" --technology "smartseq3" \
 --chemistry "SC5P-R2" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard-5pr1" --technology "smartseq3" \
 --chemistry "SC5P-R1" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

 # call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hard" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" --verbose #--as-is # "local" --localcores 1

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_hard_[IR][12]_001.fastq*
fi

#exit 0

if [[ -d input4cellranger_test-smartseq3-gelcut-small ]]; then
    rm -rf input4cellranger_test-smartseq3-gelcut-small
fi
if [[ -d test-smartseq3-gelcut-small ]]; then
    rm -rf test-smartseq3-gelcut-small
fi


# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-5pr2-hg38" --technology "smartseq3" \
 --chemistry "SC5P-R2" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-5pr1-hg38" --technology "smartseq3" \
 --chemistry "SC5P-R1" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-hg38" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "~/reference/cellranger/refdata-cellranger-GRCh38-3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-5pr2" --technology "smartseq3" \
 --chemistry "SC5P-R2" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small-5pr1" --technology "smartseq3" \
 --chemistry "SC5P-R1" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

# call on smartseq3 with files
bash launch_universc.sh --id "test-smartseq3-gelcut-small" --technology "smartseq3" \
 --chemistry "SC5P-PE" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --read1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R1.fastq" \
 --read2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_R2.fastq" \
 --index1 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I1.fastq" \
 --index2 "test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_I2.fastq" \
 --barcodefile "whitelists/Smartseq3_Fibroblasts_GelCut_small_bcs.txt" \
 --per-cell-data --jobmode "sge" #--as-is # "local" --localcores 1

# --barcodefile "whitelists/SmartSeq3_test_barcodes.txt" \

if [ -f test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_R1_001.fastq ]; then
    rename "s/_S1_L001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_S1_L001_[IR][12]_001.fastq*
    rename "s/_001//" test/shared/smartseq3-test/Smartseq3_Fibroblasts_GelCut_small_[IR][12]_001.fastq*
fi

