#!/bin/bash
#$ -N cellranger_NPO1
#$ -q bigmem.q
#$ -o cellranger_NP01.out
#$ -e cellranger_NP01_error.out
#$ -cwd

#create virtual directory of modified files
if [ ! -d cellranger ]
    then
    mkdir -p cellranger
    fi
for file in `ls *R2*.fastq | cut -d"_" -f1-3`
    do
    if [ ! -f cellranger/${file}_R1_001.fastq ]
        then
        #ln -s ../${file}_cellranger_R1_001.fastq cellranger/${file}_R1_001.fastq
        ln -s ../${file}_cellranger_R1_001.fastq cellranger/${file}_R1_001.fastq
        fi
    if [ ! -f cellranger/${file}_R2_001.fastq ]
        then
        ln -s ../${file}_R2_001.fastq cellranger/${file}_R2_001.fastq
        fi
    done

export PATH=/home/tom/local/bin/cellranger-2.1.0:$PATH

date
start=`date +%s`
cellranger count --id="test_AtRTD2" \
                 --fastqs="cellranger" \
                 --sample="NP01" \
                 --lanes="1,2" \
                 --r1-length="26" \
                 --chemistry="threeprime" \
                 --transcriptome="/home/tom/reference/cellranger/AtRTDv2/" 

#                 --noexit \
#                 --nopreflight

#                 --force-cells="2500" \
end=`date +%s`
date

runtime=$((end-start))

#edit html output
sed -i "s/test_AtRTD2/Arabidopsis \(NP01\)/g" test_AtRTD2/outs/web_summary.html
sed -i "1621s/.*/  <td><i>Arabidopsis thaliana<\/i> (NP01)<\/td>/" test_AtRTD2/outs/web_summary.html

