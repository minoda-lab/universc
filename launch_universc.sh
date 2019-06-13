#!/bin/bash
dr=''
install=false
ver_info=`paste -d "\n" <(cellranger count --version) <(echo conversion script version 0.1)`
help="Usage: bash $(basename "$0") -R1=FILE1 -R2=FILE2 -t=TECHNOLOGY [--option=OPT]
bash $(basename "$0") -v
bash $(basename "$0") -h

Convert sequencing data (FASTQ) from various platforms for compatibility with 10x Genomics

Mandatory arguments to long options are mandatory for short options too.
  -R1, --read1=FILE     Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2=FILE     Read 2 FASTQ file to pass to cellranger
  -f,  --file=NAME      Name of FASTQ files to pass to cellranger
  -t,  --technology=PLATFORM    Name of technology used to generate data (10x, nadia, icell8)
  -h,  --help     display this help and exit
  -v,  --version  output version information and exit
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
echo files: $read1 \(Read1\) and $read2 \(Read2\)

if [ "$technology" == "10x" ]; then
        echo "running cellranger for 10x"
elif [ "$technology" == "nadia" ]; then
        echo "convert barcodes for nadia"
        echo "running cellranger for nadia"
elif [ "$technology" == "icell8" ]; then
        echo "convert barcodes for iCELL8"
        echo "running cellranger for iCELL8"
fi


