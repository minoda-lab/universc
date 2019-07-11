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
bash $(basename "$0") -R1=READ1_LANE1 READ1_LANE2 -R2=READ2_LANE1 READ2_LANE2 -t=TECHNOLOGY -i=ID -r=REFERENCE [--option=OPT]
bash $(basename "$0") -f="SAMPLE_LANE" -t=TECHNOLOGY -i=ID -r=REFERENCE [--option=OPT]
bash $(basename "$0") -f="SAMPLE_LANE1" "SAMPLE_LANE2" -t=TECHNOLOGY -i=ID -r=REFERENCE [--option=OPT]
bash $(basename "$0") -v
bash $(basename "$0") -h
bash $(basename "$0") -t=TECHNOLOGY --setup

Convert sequencing data (FASTQ) from various platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
  -R1, --read1=FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2=FILE             Read 2 FASTQ file to pass to cellranger
  -f,  --file=NAME              Name of FASTQ files to pass to cellranger (prefix before R1 or R2)
  -a,  --directory=DIR          Path of directory containing all R1 and R2 files
  -t,  --technology=PLATFORM    Name of technology used to generate data (10x, nadia, icell8)
  -d,  --description=TEXT       Sample description to embed in output files.
  -c,  --chemistry=CHEM         Assay configuration, autodetection is not possible for converted files:  'SC3Pv2', 'SC5P-PE', or 'SC5P-R2' 
  -l,  --lanes=NUMS             Comma-separated lane numbers
  -n,  --force-cells=NUM        Force pipeline to use this number of cells, bypassing the cell detection algorithm.
  -i,  --id=ID                  A unique run id, used to name output folder
  -j,  --jobmode=MODE           Job manager to use. Valid options: local (default), sge, lsf, or a .template file
  -r,  --reference=DIR          Path of directory containing 10x-compatible reference.
  -h,  --help                   Display this help and exit
  -v,  --version                Output version information and exit
  -s,  --setup                  Set up whitelists for compatibility with new technology
  -k,  --keep                   Keep converted files
  -x,  --skip                   Skip conversion and run cellranger
       --verbose                Print additional outputs for debugging

For each fastq file, follow the following naming convention:
  <SampleName>_<SampleNumber>_<LaneNumber>_<ReadNumber>_001.fastq
  e.g. EXAMPLE_S1_L001_R1_001.fastq
       Example_S4_L002_R2_001.fastq.gz

Files will be renamed if they do not follow this format. File extension will be detected automatically.
"

if [[ -z $@ ]]; then
    echo "$help"
    exit 1
fi

#reading in options
read1=()
read2=()
keep=false
convert=true
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
              if [[ "$technology" != "10x" ]] && [[ "$technology" != "nadia" ]] && [[ "$technology" != "icell8" ]]; then
                echo "Error: Technology needs to be 10x, nadia, or icell8"
                exit 1
              fi
              shift
              skip=true
            fi
            ;;
        -d|--description)
            shift
            if [[ "$1" != "" ]]; then
                description="${1/%\//}"
                shift
                skip=true
            fi
            ;;
        -c|--chemistry)
            shift
            if [[ "$1" != "" ]]; then
                chemistry="${1/%\//}"
                shift
                skip=true
            fi
            ;;
        -l|--lanes)
            shift
            if [[ "$1" != "" ]]; then
                lanes="${1/%\//}"
                shift
                skip=true
            fi
            ;;
        -n|--force-cells)
            shift
            if [[ "$1" != "" ]]; then
                ncells="${1/%\//}"
                shift
                skip=true
            fi
            ;;
        -j|--jobmode)
            shift
            if [[ "$1" != "" ]]; then
                jobmode="${1/%\//}"
                shift
                skip=true
            fi
            ;;
        -i|--id)
            shift
            if [[ "$1" != "" ]]; then
                id="${1/%\//}"
                shift
                skip=true
            else
                echo "Error: identifier required"
                exit 1
            fi
            ;;
        -r|--reference)
            shift
            if [[ "$1" != "" ]]; then
                reference="${1/%\//}"
                shift
                skip=true
            else
                echo "Error: reference transcriptome generated by cellranger mkfastq required"
                exit 1
            fi
            ;;
        -f|--file)
            shift
            if [[ "$1" != "" ]]
                then
               arg=$1
               namesarray=($arg)
               for i in ${!namesarray[@]}
                    do
                    name=${namesarray[$i]}
                    echo "file: $name"
                    read1+=("${name}_R1_001")
                    read2+=("${name}_R2_001")
                    shift
                    arg=$1
                done
                skip=true
            else
                echo "Error: File input missing --file or --read1"
                exit 1
            fi
            ;;
        -a|--directory)
            shift
            if [[ "$1" != "" ]]; then
                dir=$(echo $1/*)
                for fq in $dir; do
                    if [[ $fq == *"_R1_"* ]]; then
                        read1+=("$fq")
                    elif [[ $fq == *"_R2_"* ]]; then
                        read2+=("$fq")
                    fi
                done
                shift
                skip = true
            else
                echo "Error: File input missing --directory must contain fastq files"
                exit 1
            fi
            ;;
        -R1|--read1)
            shift
            if [[ "$1" != "" ]]
                then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]
                    do
                    echo "file: $arg"
                    read1+=("${1/%\//}")
                    shift
                    arg=$1
                done
                skip=true
            else
                if [[ -z $read1 ]]; then
                    echo "Error: File missing for --read1"
                    exit 1
                fi
            fi
            ;;
        -R2|--read2)
            shift
            if [[ "$1" != "" ]]
                then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]
                    do
                    echo "file: $arg"
                    read2+=("${1/%\//}")
                    shift
                    arg=$1
                done
                skip=true
            else
                if [[ -z $read2 ]]; then
                    echo "Error: File missing for --read2"
                    exit 1
                fi
            fi
            ;;
        -s|--setup)
            echo " configure whitelist for ${cellrangerpass}..."
            setup=true
            skip=false
            shift
            ;;
        -k|--keep)
            echo " converted files will be kept if already exists"
            keep=true
            skip=false
            shift
            ;;
        -x|--skip)
            echo " running ${cellrangerpass} on existing converted files"
            keep=true
            convert=false
            skip=false
            shift
            ;;
           --verbose)
            echo " debugging mode activated"
            verbose=true
            skip=false
            ;;
        -*)
            echo "Error: Invalid option: $1"
            shift
            exit 1
            ;;
    esac
done
if [[ $verbose == "true" ]]
    then
    echo " checking options ..."
fi

#check technology input
if [[ -z $technology ]]; then
    echo "Error: option -t is required"
    exit 1
fi
#check matches expected inputs
if [[ "$technology" != "10x" ]] && [[ "$technology" != "nadia" ]] && [[ "$technology" != "icell8" ]]
    then
    echo "Error: option -t needs to be 10x, nadia, or icell8"
    exit 1
fi
#warn about barcodes for converted scripts
if [[ "$technology" == "nadia" ]] || [[ "$technology" == "icell8" ]]
    then
    echo "Warning: converting whitelist for compatibility, valid barcodes cannot be detected accurately with this technology."
fi


if [[ -z $setup ]]; then
    setup=false
    echo " skipping setup: checking if not required..."
fi

#run setup if called
if [[ $setup == "true" ]]
    then
    if [[ -z $technology ]]; then
        echo "Error: option -t is required"
        exit 1
    fi
    DIR=`which cellranger`
    VERSION=`cellranger count --version | head -n 2 | tail -n 1 | cut -d"(" -f2 | cut -d")" -f1`
    if [[ ! -w ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes ]]
        then
        echo "Error: conversion can only be run on a local install of cellranger"
        echo "Running cellranger installed at $DIR"
        echo "Install cellranger in a directory with write permissions such as /home/`whoami`/local"
        echo "Cellranger must be exported to the PATH"
        echo "The following versions of cellranger are found:"
        where cellranger
    fi
    cd ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes
    if [ -f  3M-february-2018.txt ]
        then
            gzip -f 3M-february-2018.txt
    fi
    echo "update barcodes in ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes \n for cellranger version $VERSION installed in $DIR"
    if [[ $technology == "10x" ]]
        then
        #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
        if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
            then
            echo "restore 10x barcodes"
            cp 737K-august-2016.txt.backup 737K-august-2016.txt
        fi
        echo "whitelist converted for 10x compatibility with version 2 kit"
        #create version 3 files if version 3 whitelist available
        if [ -f 3M-february-2018.txt.gz ]
            then
            #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
            if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
                then
                echo "restore 10x barcodes"
                cp 3M-february-2018.txt.backup.gz 3M-february-2018.txt.gz
            fi
            echo "whitelist converted for 10x compatibility with version 3 kit"
        fi
        #if cellranger version is 3 or greater, restore assert functions
        if [[ `printf '%s\n' '$VERSION 3.0.0' | sort -V | head -n 1` != $VERSION ]]
            then
            #restore functions if other technology run previously
            last=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called`
            if [[ $last != $technology ]]
            then
                sed -i "s/#assert barcode_idx >= prev_barcode_idx/assert barcode_idx >= prev_barcode_idx/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/assert np.array_equal(in_mc.get_barcodes(), barcodes)/g"  ${DIR}-cs/${VERSION}/lib/python/cellranger/molecule_counter.py
            fi
            echo "$DIR restored for $technology"
        else
            echo "$DIR ready for $technology"
        fi
    elif [[ $technology == "nadia" ]]
        then
        #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
        if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
            then
            echo "restore 10x barcodes"
            cp 737K-august-2016.txt.backup 737K-august-2016.txt
        fi
        #create a file with every possible barcode (permutation)
        if [ -f nadia_barcode.txt.gz ]
            then
            if [ ! -f nadia_barcode.txt ]
                then
                gunzip nadia_barcode.txt.gz
            fi
        fi
        if [ ! -f nadia_barcode.txt ]
            then
            echo AAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > nadia_barcode.txt
            echo "expected barcodes generated for Nadia"
        fi 
        #save original barcode file (if doesn't already exist)
        if [ ! -f  737K-august-2016.txt.backup ]
            then
            echo "backup of version 2 whitelist"
            cp 737K-august-2016.txt 737K-august-2016.txt.backup
        fi
        #combine 10x and Nadia barcodes
        cat nadia_barcode.txt 737K-august-2016.txt.backup | sort | uniq > 737K-august-2016.txt
        echo "whitelist converted for Nadia compatibility with version 2 kit"
        #create version 3 files if version 3 whitelist available
        if [ -f 3M-february-2018.txt.gz ]
            then
            #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
            if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
                then
                echo "restore 10x barcodes"
                cp 3M-february-2018.txt.backup.gz 3M-february-2018.txt.gz
            fi
            gunzip -k  3M-february-2018.txt.gz
            if [ ! -f  3M-february-2018.txt.backup.gz ]
                then
                echo "backup of version 3 whitelist"
                cp 3M-february-2018.txt 3M-february-2018.txt.backup
                gzip 3M-february-2018.txt.backup
            fi
            #combine 10x and Nadia barcodes
            if [ ! -f nadia_barcode.txt.gz ]
                then
                gzip -f nadia_barcode.txt
            fi
            zcat nadia_barcode.txt.gz 3M-february-2018.txt.backup.gz | sort | uniq > 3M-february-2018.txt
            gzip -f 3M-february-2018.txt
            echo "whitelist converted for Nadia compatibility with version 3 kit"
        fi
        #if cellranger version is 3 or greater, restore assert functions
        if [[ `printf '%s\n' '$VERSION 3.0.0' | sort -V | head -n 1` != $VERSION ]]
             then
             #disable functions if 10x technology run previously
             last=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called`
             if [[ $last != "10x" ]]
                 then
                 sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                 sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g"  ${DIR}-cs/${VERSION}/lib/python/cellranger/molecule_counter.py
             fi
             echo "$DIR restored for $technology"
         else
             echo "$DIR ready for $technology"
         fi
    elif [[ $technology == "icell8" ]]
        then
        #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
        if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
            then
            echo "restore 10x barcodes"
            cp 737K-august-2016.txt.backup 737K-august-2016.txt
        fi
        #create a file with every possible barcode (permutation)
        if [ -f iCELL8_barcode.txt.gz ]
            then
            if [ ! -f iCELL8_barcode.txt ]
                then
                gunzip iCELL8_barcode.txt.gz
            fi
        fi
        if [ ! -f iCELL8_barcode.txt ]
            then
            echo AAAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > iCELL8_barcode.txt
            echo "expected barcodes generated for iCELL8"
        fi
        #save original barcode file (if doesn't already exist)
        if [ ! -f  737K-august-2016.txt.backup ]
            then
            echo "backup of version 2 whitelist"
            cp 737K-august-2016.txt 737K-august-2016.txt.backup
        fi
        #combine 10x and Nadia barcodes
        cat iCELL8_barcode.txt 737K-august-2016.txt.backup | sort | uniq > 737K-august-2016.txt
        echo "whitelist converted for iCELL8 compatibility with version 2 kit"
        #create version 3 files if version 3 whitelist available
        if [ -f 3M-february-2018.txt.gz ]
            then
            #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
            if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
                then
                echo "restore 10x barcodes"
                cp 3M-february-2018.txt.backup.gz 3M-february-2018.txt.gz
            fi
            gunzip -k -f 3M-february-2018.txt.gz
            if [ ! -f  3M-february-2018.txt.backup.gz ]
                then
                echo "backup of version 3 whitelist"
                cp 3M-february-2018.txt 3M-february-2018.txt.backup
                gzip 3M-february-2018.txt.backup
            fi
            #combine 10x and Nadia barcodes
            if [ ! -f iCELL8_barcode.txt.gz ]
                then
                gzip -f iCELL8_barcode.txt
            fi
            zcat iCELL8_barcode.txt.gz 3M-february-2018.txt.backup.gz | sort | uniq > 3M-february-2018.txt
            gzip -f 3M-february-2018.txt
            echo "whitelist converted for iCELL8 compatibility with version 3 kit"        
        fi
        #if cellranger version is 3 or greater, restore assert functions
        if [[ `printf '%s\n' '$VERSION 3.0.0' | sort -V | head -n 1` != $VERSION ]]
            then
            #disable functions if 10x technology run previously
            last=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called`
            if [[ $last != "10x" ]]
                then
                sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g"  ${DIR}-cs/${VERSION}/lib/python/cellranger/molecule_counter.py
            fi
            echo "$DIR restored for $technology"
        else
            echo "$DIR ready for $technology"
        fi
    else
        echo "Error: technology not supported"
        cd -
        exit 1
    fi
    echo $technology > ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called
    cd -
    echo "setup complete"
else

#detect whitelist directory
DIR=`which cellranger`
VERSION=`cellranger count --version | head -n 2 | tail -n 1 | cut -d"(" -f2 | cut -d")" -f1`
#checkwhitelist and run --setup if needed
if [[ -f ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called ]]
    then
    #run again if last technology called is different from technology
    last=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called`
    if [[ $last != $technology ]]
        then
        echo " running setup on $technology on whitelist in ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes ..."
        if [[ $verbose == "true" ]]
            then
            echo " last: $last"
            echo " technology: $technology"
        fi
        bash $(basename "$0") -t $technology --setup
        setup=false
        if [ ! $convert ]
            then
            echo "  warning: technology changed to $technology since last run"
            convert=true
        fi
    fi
else    
        echo " using setup $technology from previous whitelist configuration ..."
        echo $technology > ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called
fi

#check inputs
if [[ -z $read1 ]]; then
    if [[ -z $file ]]; then
        echo "Error: option -R1 or --file is required"
        exit 1
    fi
elif [[ -z $read2 ]]; then
    if [[ -z $file ]]; then
        echo "Error: option -R2 or --file is required"
        exit 1
    fi
elif [[ -z $technology ]]; then
    echo "Error: option -t is required"
    exit 1
elif [[ -z $id ]]; then
    echo "Error: option -i is required"
    exit 1
elif [[ -z $reference ]]; then
    echo "Error: option -r is required"
    exit 1
fi

if [[ $id == *" "* ]]; then
    echo "id: \"$id\" must not contain a space"
    exit 1
fi

if [[ -z $description ]]; then
    description=$id
    echo "Warning: no description given, setting to ID value: $id"
fi
if [[ -z $jobmode ]]; then
    echo jobmode="local"
    echo " defaulting to local mode: --jobmode \"sge\" is recommended if running script with qsub"
fi
if [[ -z $chemistry ]]; then
    chemistry="SC3Pv2"
    echo "Warning: option -c not found, defaulting to SC3Pv2 (three-prime)"
fi
if [[ -z $keep ]]; then
    keep=false
fi

if [ -z "$ncells" ]; then
    ncells=""
else
    if ! [[ $ncells =~ '^[0-9]+$' ]] ; then
         echo "Warning --force-cells $ncells must be a number"
   fi
fi


#report inputs
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

#auto detect read file format (if one file each)
for i in ${!read1[@]}
do
read=${read1[$i]}
    echo " checking file format for $read1 ..."
    if [ -f $read ] && [ ! -h $read ]
        then
        echo $read
    elif [ -f ${read}.fq ] && [ ! -h $read ]
    then
        read=${read}.fq
        echo $read
    elif [ -f ${read}.fastq ] && [ ! -h $read ]
        then
        read=${read}.fastq
        echo $read
    elif [ -f ${read}.fq.gz ] && [ ! -h $read ]
        then
        gunzip -k ${read}.fq.gz
        read=${read}.fq
        echo $read
    elif [ -f ${read}.fastq.gz ] && [ ! -h $read ]
        then
        gunzip -k ${read}.fastq.gz
        read=${read}.fastq
        echo $read
    else
        echo $read not found
    fi
read1[$i]=$read
done

for i in ${!read2[@]}
do
read=${read2[$i]}
    echo " checking file format for $read2 ..."
    if [ -f $read ] && [ ! -h $read ]
        then
        echo $read 
    elif [ -f ${read}.fq ] && [ ! -h $read ]
    then
        read=${read}.fq
        echo $read 
    elif [ -f ${read}.fastq ] && [ ! -h $read ]
        then
        read=${read}.fastq
        echo $read 
    elif [ -f ${read}.fq.gz ] && [ ! -h $read ]
        then
        gunzip -k ${read}.fq.gz
        read=${read}.fq
        echo $read 
    elif [ -f ${read}.fastq.gz ] && [ ! -h $read ]
        then
        gunzip -k ${read}.fastq.gz
        read=${read}.fastq
        echo $read
    else
        echo $read not found
    fi
read2[$i]=$read
done

echo " checking file name for $read1 ..."
for i in ${!read1[@]}
do
read=${read1[$i]}
case $read in
    #check if contains lane before read
    (*_L0[0123456789][0123456789]_R[12]*) echo " $read compatible with lane";;
    (*) echo "  converting $read..."
        #rename file
        echo "   assuming 1 lane if none given"
        rename "s/_R1/_L001_R1/" $read
        #update file variable
        read=`echo $read | sed -e  "s/_R1/_L001_R1/g"`
        read1[$i]=$read
        echo "   renaming $read ...";;
esac
case $read in
    #check if contains sample before lane
    (*_S[123456789]_L0*) echo " $read compatible with sample";;
    (*) echo "  converting $read ..."
        #rename file
        j=$((${i}+1))
        rename "s/_L0/_S${j}_L0/" $read
        #update file variable
        read=`echo $read | sed -e  "s/_L0/_S${j}_L0/g"`
        read1[$i]=$read
        echo "   renaming $read ...";;
esac
case $read in
    #check if contains sample before lane
    (*_R1_001.*) echo " $read compatible with suffix";;
    (*) echo "  converting $read ..."
        #rename file
        rename "s/_R1.*\./_R1_001\./" $read
        #update file variable
        read=`echo $read | sed -e  "s/_R1.*\./_R1_001\./g"`
        read1[$i]=$read
        echo "   renaming $read";;
esac
#check if file in directory and link if not
path=$read
read=`echo $read | sed -r "s/.+\/(.+)/\1/"`
if [[ $verbose == "true" ]]
then
    echo "read: $read"
    echo "path: $path"
fi
if [ ! -f $read ] || [ -h $read ]
    then
    echo "  warning: file $read not in current directory"
    echo "   creating link to full path"
    ln -sf $path $read
    read1[$i]=$read
fi
done

echo " checking file name for $read2 ..."
for i in ${!read2[@]}
do
read=${read2[$i]}
case $read in
    #check if contains lane before read
    (*_L0[0123456789][0123456789]_R[12]*) echo " $read compatible with lane";;
    (*) echo "  converting $read..."
        #rename file
        echo "   assuming 1 lane if none given"
        rename "s/_R2/_L001_R2/" $read
        #update file variable
        read=`echo $read | sed -e  "s/_R2/_L001_R2/g"`
        read2[$i]=$read
        echo "   renaming $read ...";;
esac
case $read in
    #check if contains sample before lane
    (*_S[123456789]_L0*) echo " $read compatible with sample";;
    (*) echo "  converting $read ..."
        #rename file
        j=$((${i}+1))
        rename "s/_L0/_S${j}_L0/" $read
        #update file variable
        read=`echo $read | sed -e  "s/_L0/_S${j}_L0/g"`
        read2[$i]=$read
        echo "   renaming $read ...";;
esac
case $read in
    #check if contains sample before lane
    (*_R1_001.*) echo " $read compatible with suffix";;
    (*) echo "  converting $read ..."
        #rename file
        rename "s/_R1.*\./_R1_001\./" $read
        #update file variable
        read=`echo $read | sed -e  "s/_R1.*\./_R1_001\./g"`
        read2[$i]=$read
        echo "   renaming $read";;
esac
#check if file in directory and link if not
path=$read
read=`echo $read | sed -r "s/.+\/(.+)/\1/"`
if [[ $verbose == "true" ]]
then
   echo "read: $read"
   echo "path: $path"
fi
if [ ! -f $read ]  || [ -h $read ]
    then
    echo "  warning: file $read not in current directory"
    echo "   creating link to full path"
    ln -s -f $path $read
    read2[$i]=$read
fi
done

echo " files: ${read1[@]}  \(Read1\)"
echo "        ${read2[@]} \(Read2\)"

#checking the quality of fastq file names
SAMPLE=""
LANE=()

for fq in "${read1[@]}"; do
    name=`basename $fq | cut -f1 -d'.' | grep -o "_" | wc -l | xargs`
    sn=`basename $fq | cut -f1-$(($name-3))  -d'_'`
    ln=`basename $fq | cut -f$(($name-1))  -d'_' | sed 's/L00//'`
    LANE+=($ln)
    if [[ $name < 4 ]]; then
        echo "Error: filename $fq is not following the naming convention. (e.g. EXAMPLE_S1_L001_R1_001.fastq)";
        exit 1
    elif [[ $fq != *'.fastq'* ]] && [[ $fq != *'.fq'* ]]; then
        echo "Error: $fq does not have a .fq or .fastq extention"
        exit 1
    fi
    
    if [[ $sn != $SAMPLE ]]; then
        if [[ -z $SAMPLE ]]; then
            SAMPLE=$sn
        else
            echo "Error: some samples are labeled $SAMPLE while others are labeled $sn. cellranger can only handle files from one sample at a time."
            exit 1
        fi
    fi
done
for fq in "${read2[@]}"; do
    name=`basename $fq | cut -f1 -d'.' | grep -o "_" | wc -l | xargs`
    sn=`basename $fq | cut -f1-$(($name-3))  -d'_'`
    ln=`basename $fq | cut -f$(($name-1))  -d'_' | sed 's/L00//'`
    LANE+=($ln)
    
    if [[ $name < 4 ]]; then
        echo "Error: filename $fq is not following the naming convention. (e.g. EXAMPLE_S1_L001_R1_001.fastq)";
        exit 1
    elif [[ $fq != *'.fastq'* ]] && [[ $fq != *'.fq'* ]]; then
        echo "Error: $fq does not have a .fq or .fastq extention"
        exit 1
    fi
    
    if [[ $sn != $SAMPLE ]]; then
        if [[ -z $SAMPLE ]]; then
            SAMPLE=$sn
        else
            echo "Error: some samples are labeled $SAMPLE while others are labeled $sn. cellranger can only handle files from one sample at a time."
            exit 1
        fi
    fi
done

LANE=$(echo "${LANE[@]}" | tr ' ' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

if [[ $verbose == "true" ]]
    then
    echo "keep: $keep"
fi

#create directory of modified files
echo "    creating a folder for all cellranger input files..."
crIN="cellranger"
if [ ! $crIN ]
   then
   echo "    directory $crIN created for converted files"
else
   echo "    directory $crIN already exists"
fi
mkdir -p $crIN
if [ -d $crIN ]
    then
    for i in ${!read1[@]}
    do
        read=${read1[$i]}
        if [ -f $crIN/$read ] 
            then
            echo "warning: converted $read file already exists"
        fi
        if [[ $keep == "false" ]]
            then
            echo "file $crIN/$read removed"
            rm -rf $crIN/$read
        else
             echo "file $crIN/$read kept"
        fi
    done
    for i in ${!read2[@]}
    do
        read=${read2[$i]}
        if [ -f $crIN/$read ]
            then
            echo "warning: converted $read file already exists"
        fi
        if [[ $keep == "false" ]]
            then
            echo "file $crIN/$read removed"
            rm -rf $crIN/$read
        else 
             echo "file $crIN/$read kept"
        fi
    done
fi

crR1s=()
for fq in "${read1[@]}"; do
    to=`basename $fq`
    to="${crIN}/${to}"
    to=$(echo "$to" | sed 's/\.gz$//')
    crR1s+=($fq)
    
    if [[ $fq == *'.gz' ]]; then
        echo "    unzipping and redirecting $fq file..."
        gunzip -c $fq > $to
    elif [[ $keep == "false" ]]
        then
        echo "    redirecting $fq file..."
        cp $fq $to
    else
        if [ -f $crIN/$read ]
            then
            echo "    converted $fq file exists"
        else
            echo "    converted $fq file does not exist"
            cp $fq $to 
            convert=true
            echo "      warning: file will be converted"
        fi
    fi
done

crR2s=()
for fq in "${read2[@]}"; do
    to=`basename $fq`
    to="${crIN}/${to}"
    to=$(echo "$to" | sed 's/\.gz$//')
    crR2s+=($fq)
    
    if [[ $fq == *'.gz' ]]; then
        echo "    unzipping and redirecting $fq file..."
        gunzip -c $fq > $to
    elif [[ $keep  == "false" ]]
        then
        echo "    redirecting $fq file..."
        cp $fq $to
    else
        if [ -f $crIN/$read ]
            then
            echo "    converted $fq file exists"
        else
            echo "    converted $fq file does not exist"
            cp $fq $to 
            convert=true
            echo "      warning: file will be converted"
        fi
    fi
done

if [[ $verbose == "true" ]]
    then
    echo "convert: $convert"
    echo "${crR1s[@]}"
    echo "${crR2s[@]}"
fi

if [[ $convert == "true" ]]
    then
    if [[ "$technology" == "10x" ]]; then
        echo "    10x files accepted without conversion"
    else
        echo "    converting file from $technology format to 10x format..."
        for fq in "${read1[@]}"; do
            echo "        handling $fq"
            if [[ "$technology" == "nadia" ]]; then
                echo "        converting barcodes"
                sed -i '2~4s/^/AAAA/' ${crIN}/$fq #Add AAAA to every read
                echo "        converting quality scores"
                sed -i '4~4s/^/IIII/'  ${crIN}/$fq #Add quality scores for added bases
                echo "        converting UMI"
                sed -i '2~4s/[NATCG][NATCG][NATCG][NATCG][NATCG][NATCG]$/AA/'  ${crIN}/$fq #Replace last 6 bases with AA
                echo "        converting quality scores"
                sed -i '4~4s/......$/II/'  ${crIN}/$fq #Replace quality scores for added bases
            elif [[ "$technology" == "icell8" ]]; then
                echo "converting barcodes"
                sed -i '2~4s/^/AAAAA/'  ${crIN}/$fq #Add AAAAA to every read
                echo "converting quality scores"
                sed -i '4~4s/^/IIIII/'  ${crIN}/$fq #Add quality scores for added bases
            fi
        done
        echo "    conversion complete"
    fi
#close conversion
fi

#running cellranger
echo "    running cellranger..."
d=""
if [[ -n $description ]]; then
    d="--description=$description"
fi
n=""
if [[ -n $ncells ]]; then
    n="--force-cells=$ncells"
fi
j=""
if [[ -n $jobmode ]]; then
    j="--jobmode=$jobmode"
fi
start=`date +%s`
echo "cellranger count --id=$id \
        --fastqs=$crIN \
        --lanes=$LANE \
        --r1-length="26" \
        --chemistry=$chemistry \
        --transcriptome=$reference \
        --sample=$SAMPLE \
        $d \
        $n \
        $j
"

cellranger count --id=$id \
        --fastqs=$crIN \
        --lanes=$LANE \
        --r1-length="26" \
        --chemistry=$chemistry \
        --transcriptome=$reference \
        --sample=$SAMPLE \
        $d \
        $n \
        $j
#        --noexit
#        --nopreflight
end=`date +%s`
runtime=$((end-start))

#printing out log
log="
###Conversion tool log###
$ver_info

Original barcode format: ${technology} (then converted to 10x)

Runtime: ${runtime}s
"
echo "$log"

exit 0

fi
