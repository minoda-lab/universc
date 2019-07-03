#!/bin/bash
install=false

#cell renger version
cellrangerpass=`which cellranger`
if [[ -z $cellrangerpass ]]; then
	echo "cellranger command is not found."
	exit 1
fi
ver_info=`paste -d "\n" <(cellranger count --version) <(echo conversion script version 0.1) | head -n 3 | tail -n 2`
help="Usage: bash $(basename "$0") -R1=FILE1 -R2=FILE2 -t=TECHNOLOGY -i=ID -r=REFERENCE [--option=OPT]
bash $(basename "$0") -v
bash $(basename "$0") -h

Convert sequencing data (FASTQ) from various platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
  -R1, --read1=FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2=FILE             Read 2 FASTQ file to pass to cellranger
  -f,  --file=NAME              Name of FASTQ files to pass to cellranger
  -t,  --technology=PLATFORM    Name of technology used to generate data (10x, nadia, icell8)
  -i,  --id=ID	                A unique run id, used to name output folder
  -r,  --reference=DIR          Path of directory containing 10x-compatible reference.
  -h,  --help	                Display this help and exit
  -v,  --version                Output version information and exit

For each fastq file, follow the following naming convention:
  <SampleName>_<SampleNumber>_<LaneNumber>_<ReadNumber>_001.fastq
  e.g. EXAMPLE_S1_L001_R1_001.fastq
       Example_S4_L002_R2_001.fastq.gz
"
skip=false

for op in "$@";do
    if $skip;then skip=false;continue;fi
    case "$op" in
        -v|--version)
            echo "$ver_info"
            shift
            exit 0
            ;;
        -h|--help)
            echo "$help"
            shift
            exit 0
            ;;
        -t|--technology)
            shift
            if [[ "$1" != "" ]]; then
                technology=`echo "${1/%\//}" | tr '[:upper:]' '[:lower:]'`
                shift
                skip=true
            else
                echo "Error: File missing for --read1"
                exit 1
            fi
            ;;

        -f|--file)
            shift
            if [[ "$1" != "" ]]; then
                read1="${1/%\//}_R1_001"
                read2="${1/%\//}_R2_001"
                shift
                skip=true
            else
                echo "Error: File missing for --read1"
                exit 1
            fi
            ;;
        -R1|--read1)
            shift
            if [[ "$1" != "" ]]; then
                read1="${1/%\//}"
                shift
                skip=false
            else
                echo "Error: File missing for --read1"
                exit 1
            fi
            ;;
        -R2|--read2)
            shift
            if [[ "$1" != "" ]]; then
                read2="${1/%\//}"
                shift
                skip=true
            else
                echo "Error: File missing for --read1"
                exit 1
            fi
            ;;
        -*)
            echo "Error: Invalid option: $1"
            shift
            exit 1
            ;;
    esac
done

echo technology: $technology

if [ "$technology" == "10x" ]; then
        echo "running cellranger for 10x"
elif [ "$technology" == "nadia" ]; then
        echo "convert barcodes for nadia"
        echo "running cellranger for nadia"
elif [ "$technology" == "icell8" ]; then
        echo "convert barcodes for iCELL8"
        echo "running cellranger for iCELL8"
fi


if [ -f $read1 ]; then
        echo $read1
elif [ -f ${read1}.fq ]; then
        read1=${read1}.fq
        echo $read1
elif [ -f ${read1}.fastq ]; then
        read1=${read1}.fastq
        echo $read1
elif [ -f ${read1}.fq.gz ]; then
        gunzip -k ${read1}.fq.gz
        read1=${read1}.fq
        echo $read1
elif [ -f ${read1}.fastq.gz ]; then
        gunzip -k ${read1}.fq.gz
        read1=${read1}.fastq
        echo $read1
else
        echo $read1 not found
fi

if [ -f $read2 ]; then
        echo $read2
elif [ -f ${read2}.fq ]; then
        read2=${read2}.fq
        echo $read2
elif [ -f ${read2}.fastq ]; then
        read2=${read2}.fastq
        echo $read2
elif [ -f ${read2}.fq.gz ]; then
        gunzip -k ${read2}.fq.gz
        read2=${read2}.fq
        echo $read2
elif [ -f ${read2}.fastq.gz ]; then
        gunzip -k ${read2}.fq.gz
        read2=${read2}.fastq
        echo $read2
else
        echo $read2 not found
fi

echo files: $read1 \(Read1\) and $read2 \(Read2\)

