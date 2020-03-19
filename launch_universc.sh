#!/bin/bash

install=false

######convert version#####
convertversion="0.3.0.90008"
##########
cellrangerpath=`which cellranger` #location of cellranger
if [[ -z $cellrangerpath ]]; then
    echo "cellranger command is not found."
    exit 1
fi
cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
##########



#####locate launch_universc.sh, barcode sources, and other tools######
SOURCE="${BASH_SOURCE[0]}"
SDIR=""
RDIR=""

while [[ -h "$SOURCE" ]]; do #resolve $SOURCE until the file is no longer a symlink
    TARGET="$(readlink "$SOURCE")"
    if [[ $TARGET == /* ]]; then
        echo "SOURCE '$SOURCE' is an absolute symlink to '$TARGET'"
        SOURCE="$TARGET"
    else
        SDIR="$( dirname "$SOURCE" )"
        echo "SOURCE '$SOURCE' is a relative symlink to '$TARGET' (relative to '$SDIR')"
        SOURCE="$SDIR/$TARGET" #if $SOURCE is a relative symlink, we need to resolve it relative to the path where the symlink file was located
    fi
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
RDIR="$( dirname "$SOURCE" )"

if [[ $RDIR != $SDIR ]]; then
    echo "'$RDIR' resolves to '$SDIR'"
fi
echo "Running launch_universc.sh in '$SDIR'"

BARCODERECOVER=${SDIR}/RecoverBarcodes.pl
PERCELLSTATS=${SDIR}/ExtractBasicStats.pl
##########



#####usage statement#####
help='
Usage:
  bash '$(basename $0)' --testrun -t THECHNOLOGY
  bash '$(basename $0)' -t TECHNOLOGY --setup
  bash '$(basename $0)' -R1 FILE1 -R2 FILE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -R1 READ1_LANE1 READ1_LANE2 -R2 READ2_LANE1 READ2_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -f SAMPLE_LANE -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -f SAMPLE_LANE1 SAMPLE_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -v
  bash '$(basename $0)' -h

Convert sequencing data (FASTQ) from Nadia or iCELL8 platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
       --testrun                Initiates a test trun with the test dataset
  -R1, --read1 FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2 FILE             Read 2 FASTQ file to pass to cellranger
  -f,  --file NAME              Path and the name of FASTQ files to pass to cellranger (prefix before R1 or R2)
                                e.g. /path/to/files/Example_S1_L001
  
  -i,  --id ID                  A unique run id, used to name output folder
  -d,  --description TEXT       Sample description to embed in output files.
  -r,  --reference DIR          Path of directory containing 10x-compatible reference.
  -t,  --technology PLATFORM    Name of technology used to generate data.
                                Supported technologies:
                                  10x Genomics (16bp barcode, 10bp UMI): 10x, chromium (v2 or v3 automatically detected)
                                  CEL-Seq (8bp barcode, 4bp UMI): celseq
                                  CEL-Seq2 (6bp UMI, 6bp barcode): celseq2
                                  Drop-Seq (12bp barcode, 8bp UMI): nadia, dropseq
                                  iCell8 version 3 (11bp barcode, 14bp UMI): icell8 or custom
                                  inDrops version 1 (19bp barcode, 8bp UMI): indrops-v1, 1cellbio-v1
                                  inDrops version 2 (19bp barcode, 8bp UMI): indrops-v2, 1cellbio-v2
                                  inDrops version 3 (8bp barcodes, 6bp UMI): indrops-v3, 1cellbio-v3
                                  Quartz-Seq2 (14bp barcode, 8bp UMI): quartzseq2-384
                                  Quartz-Seq2 (15bp barcode, 8bp UMI): quartzseq2-1536
                                  Sci-Seq (8bp UMI, 10bp barcode): sciseq
                                  Smart-seq2-UMI, Smart-seq3 (11bp barcode, 8bp UMI): smartseq                                    
                                  SCRUB-Seq (6bp barcode, 10bp UMI): scrubseq
                                  SureCell (18bp barcode, 8bp UMI): surecell, biorad
                                Custom inputs are also supported by giving the name "custom" and length of barcode and UMI separated by "_"
                                  e.g. Custom (16bp barcode, 10bp UMI): custom_16_10
  -b,  --barcodefile FILE       Custom barcode list in plain text (with each line containing a barcode)
  
  -c,  --chemistry CHEM         Assay configuration, autodetection is not possible for converted files: 'SC3Pv2' (default), 'SC5P-PE', or 'SC5P-R2'
  -n,  --force-cells NUM        Force pipeline to use this number of cells, bypassing the cell detection algorithm.
  -j,  --jobmode MODE           Job manager to use. Valid options: 'local' (default), 'sge', 'lsf', or a .template file
       --localcores NUM         Set max cores the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --localmem NUM           Set max GB the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --mempercore NUM         Set max GB each job may use at one time.
                                    Only applies in cluster jobmodes.
  
  -p,  --per-cell-data          Generates a file with basic run statistics along with per-cell data 
  
       --setup                  Set up whitelists for compatibility with new technology and exit
       --as-is                  Skips the FASTQ file conversion if the file already exists
  
  -h,  --help                   Display this help and exit
  -v,  --version                Output version information and exit
       --verbose                Print additional outputs for debugging

For each fastq file, follow the naming convention below:
  <SampleName>_<SampleNumber>_<LaneNumber>_<ReadNumber>_001.fastq
  e.g. EXAMPLE_S1_L001_R1_001.fastq
       Example_S4_L002_R2_001.fastq.gz

For custom barcode and umi length, follow the format below:
  custom_<barcode>_<UMI>
  e.g. custom_16_10 (which is the same as 10x)

Files will be renamed if they do not follow this format. File extension will be detected automatically.
'

if [[ -z $@ ]]; then
    echo "$help"
    exit 1
fi
##########



#####options#####
#set options
lockfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.lock #path for .lock file
lastcallfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.last_called #path for .last_called
lastcall=`[[ -e $lastcallfile ]] &&  cat $lastcallfile || echo ""`
lastcall_b=`echo ${lastcall} | cut -f1 -d' '`
lastcall_u=`echo ${lastcall} | cut -f2 -d' '`
lastcall_p=`echo ${lastcall} | cut -f3 -d' '`
barcodedir=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes #folder within cellranger with the whitelist barcodes
barcodefile=""
crIN=input4cellranger #name of the directory with all FASTQ files given to cellranger
whitelistfile="outs/whitelist.txt" #name of the whitelist file added to the cellranger output
percellfile="outs/basic_stats.txt" #name of the file with the basic statistics of the run added to the cellranger output

#variable options
setup=false
convert=true
testrun=false
read1=()
read2=()
SAMPLE=""
LANE=()
id=""
description=""
reference=""
ncells=""
if [[ $technology == "10x" ]] || [[ $technology == "chromium" ]]; then
    #set default chemistry to auto detect 10x version 2 or 3
    chemistry="auto"
else
    #otherwise use version 2 configurations for other platforms
    chemistry="SC3Pv2"
fi
jobmode="local"
ncores=""
mem=""
percelldata=false

next=false
for op in "$@"; do
    if $next; then
        next=false;
        continue;
    fi
    case "$op" in
        --testrun)
            testrun=true
            next=false
            shift
            ;;
        -R1|--read1)
            shift
            if [[ "$1" != "" ]]; then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]; do
                    read1+=("${1/%\//}")
                    shift
                    arg=$1
                done
                next=true
            elif [[ -z $read1 ]]; then
                echo "Error: file input missing for --read1"
                exit 1
            fi
            ;;
        -R2|--read2)
            shift
            if [[ "$1" != "" ]]; then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]; do
                    read2+=("${1/%\//}")
                    shift
                    arg=$1
                done
                next=true
            elif [[ -z $read2 ]]; then
                echo "Error: file input missing for --read2"
                exit 1
            fi
            ;;
        -f|--file)
            shift
            if [[ "$1" != "" ]]; then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]; do
                    read1+=("${1}_R1_001")
                    read2+=("${1}_R2_001")
                    shift
                    arg=$1
                done
                skip=true
            else
                echo "Error: file input missing for --file"
                exit 1
            fi
            ;;
       -i|--id)
            shift
            if [[ "$1" != "" ]]; then
                id="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --id"
                exit 1
            fi
            ;;
        -d|--description)
            shift
            if [[ "$1" != "" ]]; then
                description="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --description"
                exit 1
            fi
            ;;
        -r|--reference)
            shift
            if [[ "$1" != "" ]]; then
                reference="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --reference: reference transcriptome generated by cellranger mkfastq required"
                exit 1
            fi
            ;;
        -t|--technology)
            shift
            if [[ $1 != "" ]]; then
                technology="${1/%\//}"
                technology=`echo "$technology" | tr '[:upper:]' '[:lower:]'`
                next=true
                shift
            else
                echo "Error: value missing for --technology"
                exit 1
            fi
            ;;
        -b|--barcodefile)
            shift
            if [[ "$1" != "" ]]; then
                barcodefile="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --barcodefile"
                exit 1
            fi
            ;;
        -c|--chemistry)
            shift
            if [[ "$1" != "" ]]; then
                chemistry="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --chemistry"
                exit 1
            fi
            ;;
        -n|--force-cells)
            shift
            if [[ "$1" != "" ]]; then
                ncells="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --force-cells"
                exit 1
            fi
            ;;
        -j|--jobmode)
            shift
            if [[ "$1" != "" ]]; then
                jobmode="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --jobmode"
                exit 1
            fi
            ;;
        --localcores)
            shift
            if [[ "$1" != "" ]]; then
                ncores="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --localcores"
                exit 1
            fi
            ;;
        --localmem)
            shift
            if [[ "$1" != "" ]]; then
                mem="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --localmem"
                exit 1
            fi
            ;;
        --mempercore)
            shift
            if [[ "$1" != "" ]]; then
                mem="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --mempercore"
                exit 1
            fi
            ;;
        -p|--per-cell-data)
            percelldata=true
            next=false
            shift
            ;;
        --setup)
            setup=true
            next=false
            shift
            ;;
        --as-is)
            convert=false
            next=false
            shift
            ;;
        -h|--help)
            echo "$help"
            exit 0
            ;;
        -v|--version)
            echo "launch_universc.sh version ${convertversion}"
            echo "cellranger version ${cellrangerversion}"
            exit 0
            ;;
       --verbose)
            echo "debugging mode activated"
            verbose=true
            next=false
            shift
            ;;
        -*)
            echo "Error: Invalid option: $op"
            exit 1
            ;;
    esac
done
##########



#####check if input maches expected formats#####
if [[ $verbose == "true" ]]; then
    echo "checking options ..."
fi

#check if this is a test run
if [[ $testrun == "true" ]]; then
    if [[ ${#read1[@]} -gt 0 ]] || [[ ${#read2[@]} -gt 0 ]]; then
        echo "Error: for test run, no R1 or R2 file can be selected."
        exit 1
    fi
    
    reference=${SDIR}/test/cellranger_reference/cellranger-tiny-ref/3.0.0
    
    percelldata=true
    id=test-tiny-${technology}
    description=${id}
    
    if [[ $technology == "10x" ]]; then
        gunzip -fk ${SDIR}/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L00[12]_R[12]_001.fastq.gz
        read1=("${SDIR}/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq" "${SDIR}/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq")
        read2=("${SDIR}/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_002.fastq" "${SDIR}/test/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R2_001.fastq")
    elif [[ $technology == "nadia" ]]; then
        gunzip -fk ${SDIR}/test/shared/dropseq-test/SRR1873277_S1_L001_R[12]_001.fastq.gz
        read1=("${SDIR}/test/shared/dropseq-test/SRR1873277_S1_L001_R1_001.fastq")
        read2=("${SDIR}/test/shared/dropseq-test/SRR1873277_S1_L001_R2_001.fastq")
    elif [[ $technology == "icell8" ]]; then
        gunzip -fk ${SDIR}/test/shared/mappa-test/test_FL_R[12].fastq.gz
        read1=("${SDIR}/test/shared/mappa-test/test_FL_R1.fastq")
        read2=("${SDIR}/test/shared/mappa-test/test_FL_R2.fastq")
    else
        echo "Error: for test run, option --technology must be 10x, nadia, or icell8"
	exit 1
    fi
fi

#check if cellranger is writable
if ! [[ -w "$barcodedir" ]]; then
    echo "Error: Trying to run cellranger installed at ${cellrangerpath}"
    echo "launch_universc.sh can only run with cellranger installed locally"
    echo "Install cellranger in a directory with write permissions such as /home/`whoami`/local and export to the PATH"
    echo "The following versions of cellranger are found:"
    echo " `whereis cellranger`"
    exit 1
fi

#check if convert is writable
if ! [[ -w "$SDIR" ]]; then
    echo "Error: Trying to run launch_universc.sh installed at $SDIR"
    echo "$SDIR must be writable to run launch_universc.sh"
    echo "Install launch_universc.sh in a directory with write permissions such as /home/`whoami`/local and export to the PATH"
    exit 1
fi


#aliases for technology with the same settings
if [[ "$technology" == "chromium" ]]; then
    echo "Running with technology 10x (chromium)" 
    technology="10x"
fi
if [[ "$technology" == "celseq" ]] || [[ "$technology" == "cel-seq" ]]; then
    echo "Running with CEL-Seq parameters (version 1)"
    technology="celseq"
fi
if [[ "$technology" == "celseq2" ]] || [[ "$technology" == "cel-seq2" ]]; then
    echo "Running with CEL-Seq2 parameters"
    technology="celseq2"
fi
if [[ "$technology" == "dropseq" ]] || [[ "$technology" == "drop-seq" ]]; then
    echo "Running with Nadia parameters (Drop-Seq)"
    technology="nadia"
fi
if [[ "$technology" == "icell8" ]] || [[ "$technology" == "icell-8" ]]; then
    echo "Running with iCELL8 parameters (version 3 with UMIs)"
    echo "***WARNING: iCELL8 settings should only be used for kits that have UMIs***"
    technology="icell8"
fi
if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrops-v1" ]] || [[ "$technology" == "indropv1" ]] || [[ "$technology" == "indropsv1" ]] || [[ "$technology" == "1cellbio-v1" ]] || [[ "$technology" == "1cellbiov1" ]]; then
    echo "Running with inDrop parameters (version 1)"
    technology="indrop-v1"
fi
if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrops-v2" ]] || [[ "$technology" == "indropv2" ]] || [[ "$technology" == "indropsv2" ]] || [[ "$technology" == "1cellbio-v2" ]] || [[ "$technology" == "1cellbiov2" ]]; then
    echo "Running with inDrop parameters (version 2 with reads inverted)"
    technology="indrop-v2"    
fi
if [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "indrops-v3" ]] || [[ "$technology" == "indropv3" ]] || [[ "$technology" == "indropsv3" ]] || [[ "$technology" == "1cellbio-v3" ]] || [[ "$technology" == "1cellbiov3" ]]; then
    echo "Running with inDrop parameters (version 3 with reads inverted)"
    technology="indrop-v3"
fi
if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
    #invert read1 and read2
    echo "Using barcodes on Read 2"
    tmp=$read1
    read1=$read2
    read2=$tmp
    tmp=""
fi
if [[ "$technology" == "quartz-seq2-384" ]] || [[ "$technology" == "quartzseq2-384" ]] || [[ "$technology" == "quartz-seq2-v3.1" ]] || [[ "$technology" == "quartzseq2-v3.1" ]] || [[ "$technology" == "quartzseq2v3.1" ]]; then
    echo "Running with Quartz-Seq2 v3.1 parameters 384 wells (14bp barcode)"
    technology="quartz-seq2-384"
fi

if [[ "$technology" == "quartz-seq2-1536" ]] || [[ "$technology" == "quartzseq2-1536" ]] || [[ "$technology" == "quartz-seq2-v3.2" ]] || [[ "$technology" == "quartzseq2-v3.2" ]] || [[ "$technology" == "quartzseq2v3.2" ]]; then
    echo "Running with Quartz-Seq2 v3.2 parameters 1536 wells (15bp barcode)"
    technology="quartz-seq2-1536"
fi
if [[ "$technology" == "sciseq" ]] || [[ "$technology" == "sci-seq" ]]; then
    echo "Running with Sci-Seq parameters (single-cell combinatorial indexing RNA sequencing)"
    technology="sciseq"
fi
if [[ "$technology" == "scrubseq" ]] || [[ "$technology" == "scrub-seq" ]]; then
    echo "Running with SCRUB-Seq parameters"
    technology="scrubseq"
fi
if [[ "$technology" == "smartseq" ]] || [[ "$technology" == "smart-seq" ]] || [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smart-seq2" ]] ||  [[ "$technology" == "smartseq2-umi" ]] || [[ "$technology" == "smart-seq2-umi" ]] ||  [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "smart-seq3" ]]; then
    echo "Running with Smart-Seq3 parameters (version 3 with UMIs)"
    echo "***WARNING: Smart-Seq settings should only be used for kits that have UMIs***"
    technology="smartseq"
fi
if [[ "$technology" == "surecell" ]] || [[ "$technology" == "surecellseq" ]] || [[ "$technology" == "surecell-seq" ]]|| [[ "$technology" == "bioraad" ]]; then
    echo "Running with SureCell parameters"
    technology="surecell"
fi


if [[ "$technology" == "indrop-v3" ]] \
|| [[ "$technology" == "sciseq" ]] \
|| [[ "$technology" == "smart-seq"* ]]; then
    ## see the following issues from supporting inDrops-v3
    # https://github.com/alexdobin/STAR/issues/825
    # https://github.com/BUStools/bustools/issues/4
    ## these workarounds could be implemented in the future
    echo "***WARNING: dual indexes are not supported***"
    echo "...samples should be demultiplexed into separate FASTQ files by index before running"
fi


#check if technology matches expected inputs
if [[ "$technology" != "10x" ]] \
&& [[ "$technology" != "celseq"* ]] \
&& [[ "$technology" != "nadia" ]] \
&& [[ "$technology" != "icell8" ]] \
&& [[ "$technology" != "indrop"* ]] \
&& [[ "$technology" != "quartz-seq2"* ]] \
&& [[ "$technology" != "sciseq" ]] \
&& [[ "$technology" != "scrubseq" ]] \
&& [[ "$technology" != "smart-seq"* ]]\
&& [[ "$technology" != "surecell" ]]; then
    if [[ "$technology" != "custom"* ]]; then
        echo "Error: option -t needs to be 10x, nadia, icell8, or custom_<barcode>_<UMI>"
        exit 1
    else
        custom=`echo $technology | grep -o "_" | wc -l`
        custom=$(($custom+1))
        if [[ $custom -le 2 ]]; then
           echo "Error: custom input must have at least 3 fields, e.g., customA_10_16"
        fi
        b=`echo $technology | cut -f$((${custom}-1))  -d'_'`
        u=`echo $technology | cut -f$((${custom}))  -d'_'`
        if ! [[ "$b" =~ ^[0-9]+$ ]] || ! [[ "$u" =~ ^[0-9]+$ ]]; then
            echo "Error: option -t needs to be 10x, nadia, icell8, or custom_<barcode>_<UMI>"
            exit 1
        fi
        if [[ -z $barcodefile ]]; then
            echo "Warning: when option -t is set as custom, a file with a list of barcodes can be specified with option -b."
        fi
    fi
fi

#check for presence of read1 and read2 files
if [[ $setup == "false" ]]; then
    if [[ ${#read1[@]} -eq 0 ]]; then
        echo "Error: option -R1 or --file is required"
        exit 1
    elif [[ ${#read2[@]} -eq 0 ]]; then
        echo "Error: option -R2 or --file is required"
        exit 1
    fi
fi

#check read1 and read2 files for their extensions
##allows incomplete file names and processing compressed files
for i in {1..2}; do
    readkey=R$i
    list=""
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    fi
    
    for j in ${!list[@]}; do
        read=${list[$j]}
        if [[ $verbose == "true" ]]; then
            echo "checking file format for $read ..."
        fi
        if [[ -f $read ]] && [[ -h $read ]]; then
            if [[ $read == *"gz" ]]; then
                gunzip -f -k $read
                #update file variable
                read=`echo $read | sed -e "s/\.gz//g"`
            fi
            if [[ $read != *"fastq" ]] && [[ $read != *"fq" ]]; then
                echo "Error: file $read needs a .fq or .fastq extention."
                exit 1
            fi
        elif [[ -f $read ]]; then
            if [[ $read == *"gz" ]]; then
                gunzip -f -k $read
                #update file variable
                read=`echo $read | sed -e "s/\.gz//g"`
            fi
            if [[ $read != *"fastq" ]] && [[ $read != *"fq" ]]; then
                echo "Error: file $read needs a .fq or .fastq extention."
                exit 1
            fi
        #allow detection of file extension (needed for --file input)
        elif [[ -f ${read}.fq ]] || [[ -h ${read}.fq ]]; then
            read=${read}.fq
        elif [[ -f ${read}.fastq ]] || [[ -h ${read}.fastq ]]; then
            read=${read}.fastq
        elif [[ -f ${read}.fq.gz ]] || [[ -h ${read}.fq.gz ]]; then
            gunzip -f -k ${read}.fq.gz
            read=${read}.fq
        elif [[ -f ${read}.fastq.gz ]] || [[ -h ${read}.fastq.gz ]]; then
            gunzip -f -k ${read}.fastq.gz
            read=${read}.fastq
        else
            echo "Error: $read not found"
            exit 1
        fi
        
        if [[ $verbose == "true" ]]; then
             echo "  $read"
        fi
        
        list[$j]=$read
    done
    
    if [[ $readkey == "R1" ]]; then
        read1=("${list[@]}")
    elif [[ $readkey == "R2" ]]; then
        read2=("${list[@]}")
    fi
done

#renaming read1 and read2 files if not compatible with the convention
for i in {1..2}; do
    readkey=R$i
    list=""
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    fi
    
    for j in ${!list[@]}; do
        read=${list[$j]}
        if [[ $verbose == "true" ]]; then
            echo " checking file name for $read ..."
        fi
        
        if [[ -h $read ]]; then
            path=`readlink -f $read`
            if [[ $verbose == "true" ]]; then
                echo " ***Warning: file $read not in current directory. Path to the file captured instead.***"
                echo "  (file) $read"
                echo "  (path) $path"
            fi
            read=${path}
        fi
        case $read in
            #check if contains lane before read
            *_L0[0123456789][0123456789]_$readkey*)
                if [[ $verbose == "true" ]]; then
                    echo "  $read compatible with lane"
                fi
            ;;
            *) 
                #rename file
                if [[ $verbose == "true" ]]; then
                    echo "***Warning: file $read does not have lane value in its name. Lane 1 is assumed.***"
                echo "  renaming $read ..."
                fi
                rename "s/_$readkey/_L001_$readkey/" $read
                #update file variable
                read=`echo $read | sed -e "s/_${readkey}/_L001_${readkey}/g"`
                list[$j]=$read
            ;;
        esac
        case $read in
            #check if contains sample before lane
            *_S[0123456789]_L0*)
                if [[ $verbose == "true" ]]; then
                    echo "  $read compatible with sample"
                fi
            ;;
            *)
                #rename file
                if [[ $verbose == "true" ]]; then
                    echo "***Warning: file $read does not have sample value in its name. Sample $k is assumed.***"
                    echo "  renaming $read ..."
                fi
                k=$((${j}+1))
                rename "s/_L0/_S${k}_L0/" $read
                #update file variable
                read=`echo $read | sed -e "s/_L0/_S${k}_L0/g"`
                list[$j]=$read
            ;;
        esac
        case $read in
            #check if contains sample before lane
            *_${readkey}_001.*)
                if [[ $verbose == "true" ]]; then
                    echo "  $read compatible with suffix"
                fi
            ;;
            *)
                #rename file
                if [[ $verbose == "true" ]]; then
                    echo "***Warning: file $read does not have suffix in its name. Suffix 001 is given.***"
                    echo "  renaming $read ..."
                fi
                rename "s/_${readkey}.*\./_${readkey}_001\./" $read
                #update file variable
                read=`echo $read | sed -e "s/_${readkey}.*\./_${readkey}_001\./g"`
                list[$j]=$read
            ;;
        esac
    done
    
    if [[ $readkey == "R1" ]]; then
        read1=("${list[@]}")
    elif [[ $readkey == "R2" ]]; then
        read2=("${list[@]}")
    fi
done

#checking the quality of fastq file names
read12=("${read1[@]}" "${read2[@]}")
for fq in "${read12[@]}"; do
    name=`basename $fq`
    name=${name%.*}
    fields=`echo ${name} | grep -o "_" | wc -l`
    fields=$(($fields+1))
    sn=`echo ${name} | cut -f1-$((${fields}-4))  -d'_'`
    lane=`echo ${name} | cut -f$((${fields}-2)) -d'_' | sed 's/L00//'`
    LANE+=($lane)
    if [[ ${fields} -le 4 ]]; then
        echo "Error: filename $fq is not following the naming convention. (e.g. EXAMPLE_S1_L001_R1_001.fastq)";
        exit 1
    elif [[ $fq != *'.fastq'* ]] && [[ $fq != *'.fq'* ]]; then
        echo "Error: $fq does not have a .fq or .fastq extention."
        exit 1
    elif [[ ${sn} =~ "." ]]; then
        echo "Error: $fq has a period \".\" within its sample name. Remove it to run cellranger."
	exit 1
    fi
    
    if [[ ${sn} != $SAMPLE ]]; then
        if [[ -z $SAMPLE ]]; then
            SAMPLE=${sn}
        else
            echo "Error: some samples are labeled $SAMPLE while others are labeled $sn. cellranger can only handle files from one sample at a time."
            exit 1
        fi
    fi
done
LANE=$(echo "${LANE[@]}" | tr ' ' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')

#select the input barcode file
if [[ ! -z "$barcodefile" ]]; then
    if [[ ! -f $barcodefile ]]; then
        echo "Error: File selected for option --barcodefile does not exist"
	exit 1
    else
        barcodefile=`readlink -f $barcodefile`
        #all barcodes upper case
        sed -i 's/.*/\U&/g' $barcodefile
    fi
else
    if [[ "$technology" == "10x" ]]; then
        barcodefile="default"
    elif [[ "$technology" == "nadia" ]]; then
        barcodefile=${SDIR}/nadia_barcode.txt
        if [[ ! -f ${barcodefile} ]]; then
            #create a nadia barcode file
            echo "No barcodes whitelists available for Drop-Seq or Nadia: all possible barcodes accepted (valid barcodes will be 100% as a result)"
            echo {A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' | sort | uniq > ${barcodefile}
        fi
    elif [[ "$technology" == "icell8" ]]; then
        barcodefile=${SDIR}/iCell8_barcode.txt
    elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
        # use bustools whitelist for inDrops-v2 with adapters removed https://github.com/BUStools/bustools/issues/4 
        barcodefile=${SDIR}/inDrops_barcodes.txt
    elif [[ "$technology" == "indrop-v3" ]];
        # inDrops-v3 whitelist is a combination of v2 whitelists https://github.com/indrops/indrops/issues/32
        ## version 2 whitelist will be used until dual indexing (i7) is supported for read I2: https://github.com/alexdobin/STAR/issues/825
        barcodefile=${SDIR}/inDrops-v2_barcodes.txt 
    elif [[ "$technology" == "quartz-seq2-384" ]]; then
        barcodefile=${SDIR}/Quartz-Seq2-384_barcode.txt
    elif [[ "$technology" == "quartz-seq2-1536" ]]; then
        barcodefile=${SDIR}/Quartz-Seq2-1536_barcode.txt
    elif [[ "$technology" == "smartseq" ]]; then
        echo "***WARNING: barcodes not available for Smart-Seq 2 or 3, using iCELL8 whitelist (version 3)***"
        echo "...valid barcodes may be an overestimate"
        barcodefile=${SDIR}/iCell8_barcode.txt
    elif [[ "$technology" == "custom"* ]] || [[ "$technology" == "celseq"* ]] ||  [[ "$technology" == "scrubseq" ]] || [[ "$technology" == "sciseq" ]] || [[ "$technology" == "surecell" ]]; then
        if [[ "$technology" == "celseq" ]]; then
            customname="celseq"
            minlength=8
        elif [[ "$technology" == "celseq2" ]]; then
            customname="celseq"
            minlength=6
        elif [[ "$technology" == "sciseq" ]]; then
            customname="sciseq"
            minlength=10
        elif [[ "$technology" == "scrubseq" ]]; then
            customname="scrubseq"
            minlength=6
        elif [[ "$technology" == "surecell" ]]; then
            customname="surecell"
            barcodelength=18
            minlength=$(( $barcodelength < 16 ? $barcodelength : 16 ))
        elif [[ "$technology" == "custom"* ]]; then
            custom=`echo $technology | grep -o "_" | wc -l`
            custom=$(($custom+1))
            customname=`echo $technology | cut -f1-$((${custom}-2))  -d'_'`
            barcodelength=`echo $technology | cut -f$((${custom}-1))  -d'_'`
            #check whether barcodes exceed 16bp (and reuse 16bp whitelist for all greater than 16bp)
            if [[  $barcodelength -ge 16 ]]; then
                echo "Barcode length ($barcodelength) of 16 or more:"
                echo "    ...using barcode whitelist of 16bp"
            fi
            minlength=$(( $barcodelength < 16 ? $barcodelength : 16 ))
        fi
        # compute custom barcodes if barcode length is different
        barcodefile=${SDIR}/${customname}_${minlength}_barcode.txt
        if [[ ! -f ${barcodefile} ]]; then
            echo "No barcodes whitelists available for ${customname}: all possible barcodes accepted (valid barcodes will be 100% as a result)"
            echo "***Warning: giving a barcode whitelist --barcodefile is recommended where available.***"
            if [[ -f ${SDIR}/*${barcodelength}_barcode.txt ]]; then
                pregeneratedfile=`ls ${SDIR}/*${barcodelength}_barcode.txt | awk '{print $1;}' | head -n 1`
                echo "$pregeneratedfile for barcode(${barcodelength}) generated already"
                ln -s $pregeneratedfile $barcodefile
                echo "    ...using this as $barcodefile"
            else
                echo "generating $barcodefile of barcode length $minlength"
                # generating permutations of ATCG of barcode length (non-standard evaluation required to run in script)
                echo $(eval echo $(for ii in $(eval echo {1..${minlength}}); do echo "{A,T,C,G}"; done |  tr "\n" " " | sed "s/ //g" | xargs -I {} echo {})) | sed 's/ /\n/g' | sort | uniq > ${barcodefile}
            fi
        fi 
    fi
fi

#check if reference is present
if [[ -z $reference ]]; then
    if [[ $setup == "false" ]] || [[ ${#read1[@]} -ne 0 ]] || [[ ${#read2[@]} -ne 0 ]]; then
        echo "Error: option --reference is required";
        exit 1
    fi
fi

#check if ncells is an integer
int='^[0-9]+$'
if [[ -z "$ncells" ]]; then
    ncells=""
elif ! [[ $ncells =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --force-cells must be an integer"
    exit 1
fi

#check if ncores is an integer
int='^[0-9]+$'
if [[ -z "$ncores" ]]; then
    ncores=""
elif ! [[ $ncores =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --localcores must be an integer"
    exit 1
fi

#check if mem is a number
int='^[0-9]+([.][0-9]+)?$'
if [[ -z "$mem" ]]; then
    mem=""
elif ! [[ $mem =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --localmem or --mempercore must be a number (of GB)"
    exit 1
fi

#check if chemistry matches expected input
#allow "auto" for 10x
if [[ "$technology" != "10x" ]]; then
    #use SC3Pv3 (umi length 12) 
    if [[ $umilength -ge 11 ]] && [[ "$chemistry" == "SC3Pv1" ]] && [[ "$chemistry" == "SC3Pv2" ]]; then
        echo "Using 10x version 3 chemistry to support longer UMIs"
        chemistry="SC3Pv3"
    else
    #use SC3Pv2 (umi length 10)
         echo "Using 10x version 2 chemistry to support UMIs"
        chemistry="SC3Pv2"
    fi
    if [[ "$chemistry" != "SC3Pv1" ]] && [[ "$chemistry" != "SC3Pv2" ]] && [[ "$chemistry" != "SC3Pv3" ]] && [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R2" ]]; then
        echo "Error: option --chemistry must be SC3Pv3, SC3Pv2, SC5P-PE , or SC5P-R2"
       exit 1
    fi
fi

#checking if jobmode matches expected input
if [[ "$jobmode" != "local" ]] && [[ "$jobmode" != "sge" ]] && [[ "$jobmode" != "lsf" ]] && [[ "$jobmode" != *"template" ]]; then
    echo "Error: option --jobmode must be local, sge, lsf, or a .template file"
    exit 1
fi

#check if ID is present
if [[ -z $id ]]; then
    if [[ ${#read1[@]} -ne 0 ]] || [[ ${#read2[@]} -ne 0 ]]; then
        echo "Error: option --id is required"
        exit 1
    fi
fi
crIN=${crIN}_${id}

#checking if crIN exists
if [[ ! -d $crIN ]]; then
    convert=true
    echo "***Warning: convertion was turned on because directory $crIN was not found***"
fi
##########



#####Get barcode/UMI length#####  
barcode_default=16
if [[ "$chemistry" == *"v2" ]]; then
    umi_default=10
elif [[ "$chemistry" == *"v3" ]]; then
    umi_default=12
fi

totallength=`echo $((${barcode_default}+${umi_default}))`

#barcode and umi lengths given by options
barcodelength=""
umilength=""
if [[ "$technology" == "10x" ]]; then
    barcodelength=16
    umilength=10
    umi_default=10
elif [[ "$technology" == "celseq" ]]; then
    barcodelength=8
    umilength=4
elif [[ "$technology" == "celseq2" ]]; then
    barcodelength=6
    umilength=6
elif [[ "$technology" == "nadia" ]]; then
    barcodelength=12
    umilength=8
elif [[ "$technology" == "icell8" ]]; then
    barcodelength=11
    umilength=14
elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
    barcodelength=19 
    umilength=6
elif [[ "$technology" == "indrop-v3" ]]; then
    barcodelength=8
    umilength=6
elif [[ "$technology" == "quartz-seq2-384" ]]; then
    barcodelength=14
    umilength=8
elif [[ "$technology" == "quartz-seq2-1536" ]]; then
    barcodelength=15
    umilength=8 
elif [[ "$technology" == "sciseq" ]]; then
    barcodelength=10
    umilength=8
elif [[ "$technology" == "scrubseq" ]]; then
    barcodelength=6 
    umilength=10
elif [[ "$technology" == "smartseq" ]]; then
    barcodelength=11
    umilength=8
elif [[ "$technology" == "surecell" ]]; then
    barcodelength=18
    umilength=8
else
    custom=`echo $technology | grep -o "_" | wc -l`
    custom=$(($custom+1))
    barcodelength=`echo $technology | cut -f$((${custom}-1))  -d'_'`
    umilength=`echo $technology | cut -f$((${custom}))  -d'_'`
fi


#adjustment lengths
barcodeadjust=`echo $(($barcodelength-$barcode_default))`
umiadjust=`echo $(($umilength-$umi_default))`
##########



#####check if UniverSC is running already#####
#set up .lock file
if [[ ! -f $lockfile ]]; then
    echo "creating .lock file"
    echo 0 > $lockfile
    lock=`cat $lockfile`
else
    #check if jobs are running (check value in .lock file)
    echo "checking .lock file"
    lock=`cat $lockfile`
    
    if [[ $lock -le 0 ]]; then
        echo " call accepted: no other cellranger jobs running"
        lock=1
        if [[ $setup == "false" ]]; then 
                    echo $lock > $lockfile
        fi
    else
        if [[ -f $lastcallfile ]]; then
	    echo " total of $lock cellranger ${cellrangerversion} jobs are already running in ${cellrangerpath} with barcode length (${lastcall_b}), UMI length (${lastcall_u}), and whitelist barcodes (${lastcall_p})"
            
	    #check if a custom barcode is used for a run (which cannot be run in parallel)
            if [[ ${barcodelength} == ${lastcall_b} ]] && [[ ${umilength} == ${lastcall_u} ]] && [[ ${barcodefile} == ${lastcall_p} ]]; then
                echo " call accepted: no conflict detected with other jobs currently running"
                #add current job to lock
                lock=$(($lock+1))
                if [[ $setup == "false" ]]; then 
                    echo $lock > $lockfile
                fi
            else
                echo "Error: conflict between technology selected for the new job and other jobs currently running"
                echo "make sure that the barcode length, UMI length, and the whitelist barcodes are the same as the other jobs currently running"
                echo "if confident that no other jobs are running and still get this error, remove $lockfile and try again"
                if [[ $verbose ]]; then
                    echo "Submitted configuration with barcode length (${barcodelength}), UMI length (${umilength}), and whitelist barcodes (${barcodefile})"
                fi
                exit 1
            fi
        else
            echo "Error: $lastcallfile not found"
        fi
    fi
fi
##########



####report inputs#####
echo ""
echo "#####Input information#####"
echo "SETUP and exit: $setup"
if [[ $setup == "true" ]]; then
    echo "***Warning: launch_universc.sh will exit once whitelist is converted***"
fi
echo "FORMAT: $technology"
if [[ $technology == "nadia" ]]; then
    echo "***Warning: whitelist is converted for compatibility with $technology, valid barcodes cannot be detected accurately with this technology***"
fi
echo "BARCODES: ${barcodefile}"
if [[ ${#read1[@]} -eq 0 ]] && [[ ${#read1[@]} -eq 0 ]]; then
    echo "***Warning: no FASTQ files were selected, launch_universc.sh will exit after setting up the whitelist***"
fi
if ! [[ ${#read1[@]} -eq 0 ]]; then
    echo "INPUT(R1):"
    for i in ${!read1[@]}; do
        echo " ${read1[$i]}"
    done
fi
if ! [[ ${#read2[@]} -eq 0 ]]; then
    echo "INPUT(R2):"
    for i in ${!read2[@]}; do
        echo " ${read2[$i]}"
    done
fi
echo "SAMPLE: $SAMPLE"
echo "LANE: $LANE"
echo "ID: $id"
if [[ -z $description ]]; then
    description=$id
    echo "DESCRIPTION: $description"
    echo "***Warning: no description given, setting to ID value***"
else
    echo "DESCRIPTION: $description"
fi
echo "REFERENCE: $reference"
if [[ -z $ncells ]]; then
    echo "NCELLS: $ncells(no cell number given)"
else
    echo "NCELLS: $ncells"
fi
echo "CHEMISTRY: $chemistry"
echo "JOBMODE: $jobmode"
if [[ "$jobmode" == "local" ]]; then
    echo "***Warning: --jobmode \"sge\" is recommended if running script with qsub***"
fi
echo "CONVERSTION: $convert"
if [[ $convert == "false" ]]; then
    echo "***Warning: adjustment for barcode and UMI length was skipped***"
fi
echo "##########"
echo ""
##########



####setup whitelist#####
if [[ $lock -eq 0 ]]; then
    echo "setup begin"
    echo "updating barcodes in $barcodedir for cellranger version ${cellrangerversion} installed in ${cellrangerpath} ..."
    
    cd $barcodedir
    
    #restore assert functions if cellranger version is 3 or greater
    echo " restoring cellranger"
    if [[ `printf '%s\n' '${cellrangerversion} 3.0.0' | sort -V | head -n 1` != ${cellrangerversion} ]]; then
        if [[ $technology == "10x" ]] && [[ -z $barcodefile ]]; then
            sed -i "s/#if gem_group == prev_gem_group/if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#assert barcode_idx >= prev_barcode_idx/assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        elif [[ $lastcall == "10x" ]] || [[ ! -f $lastcallfile ]]; then
            sed -i "s/if gem_group == prev_gem_group/#if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        fi
        echo " ${cellrangerpath} set for $technology"
    fi
    
    #whitelist file name
    v2=737K-august-2016.txt
    v3=3M-february-2018.txt

    #generate backup for the default 10x whitelist
    if [[ ! -f 737K-august-2016.txt.backup ]] || [[ ! -f 3M-february-2018.txt.backup.gz ]]; then
        echo " generating backups for default 10x whitelist"
        cp -f ${v2} ${v2}.backup
       	cp -f ${v3}.gz ${v3}.backup.gz
        echo " backup generated"
    fi
    
    #convert whitelist to the apropriate barcode
    echo " converting whitelist"
    if [[ ${barcodefile} == "default" ]]; then
        #for version 2
        cp ${v2}.backup ${v2}
        #for version 3
        cp ${v3}.backup.gz ${v3}.gz
    else
        #for version 2
        cat ${barcodefile} > ${v2}
        if [[ $barcodeadjust -gt 0 ]]; then
            sed -i "s/^.{${barcodeadjust}}//" ${v2} #Trim the first n characters from the beginning of the sequence and quality
        elif [[ 0 -gt $barcodeadjust ]]; then
            As=`printf '%0.sA' $(seq 1 $(($barcodeadjust * -1)))`
            sed -i "s/^/$As/" ${v2} #Trim the first n characters from the beginning of the quality
        fi
        #for version 3
        cat ${v2} > ${v3}
        gzip -f ${v3}
        rm translation/${v3}.gz
        ln -s ${v3}.gz translation/${v3}.gz
    fi
    echo " whitelist converted"
    
    #change last call file
    echo "${barcodelength} ${umilength} ${barcodefile}" > $lastcallfile
    
    cd - > /dev/null
    
    echo "setup complete"
    if [[ $setup == "true" ]]; then
        exit 0
    fi
fi
#########



#####exiting when setup is all that is requested#####
if [[ $setup == "ture" ]]; then
    lock=`cat $lockfile`
    lock=$(($lock-1))
    echo $lock > $lockfile
    echo " setup complete. exiting launch_universc.sh"
    exit 0
fi
##########



#####create directory with files fed to cellranger#####
echo "creating a folder for all cellranger input files ..."
convFiles=()

if [[ ! -d $crIN ]]; then
    echo " directory $crIN created for converted files"
    mkdir $crIN
else
    echo " directory $crIN already exists"
fi

if [[ $convert == "true" ]]; then
    echo "moving file to new location"
fi

crR1s=()
for fq in "${read1[@]}"; do
    to=`basename $fq`
    to="${crIN}/${to}"

    #invert read1 and read2
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        #where converted "read1" is R2 in source files (corrected names for cellranger)
        to=`echo $to | sed -e "s/_R2_/_R1_/g"
    fi

    crR1s+=($to)
    
    echo " handling $fq ..."
    if [[ ! -f $to ]] || [[ $convert == "true" ]]; then
        cp -f $fq $to
    fi
    if [[ $convert == "true" ]]; then
        convFiles+=($to)
    fi
done

crR2s=()
for fq in "${read2[@]}"; do
    to=`basename $fq`
    to="${crIN}/${to}"
    to=$(echo "$to" | sed 's/\.gz$//')

    #invert read1 and read2
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        #where converted "read2" is R1 in source files
        to=`echo $to | sed -e "s/_R1_/_R2_/g"
    fi

    crR2s+=($to)
    
    echo " handling $fq ..."
    if [[ ! -f $to ]] || [[ $convert == "true" ]]; then
        cp -f $fq $to
    fi
done
##########



#####convert file format#####
echo "converting input files to confer cellranger format ..."
if [[ $convert == "false" ]]; then
    echo " input file format conversion skipped"
else
    echo " adjustment parameters:"
    echo "  barcodes: ${barcodeadjust}bps at its head"
    echo "  UMIs: ${umiadjust}bps at its tail" 
    
    #for CEL-Seq2 swap barcode and UMI
    if [[ "$technology" == "sciseq" ]]; then
        for convFile in "${convFiles[@]}"; do
            #swap UMI and barcode
            sed -E '2~2s/(.{6})(.{6})/\2\1/' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi

    #remove adapter from inDrops
    if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]];
        for convFile in "${convFiles[@]}"; do
            #remove adapter if present
            sed -E '
                /^(.{8})GAGTGATTGCTTGTGACGCCTT(.{8})/ {
                s/^(.{8})GAGTGATTGCTTGTGACGCCTT(.{8})/\1\2/g
                n
                n
                s/^(.{8}).{22}(.{8})/\1\2/g
                }' $convFile |
            #remove linker between barcode and UMI
            sed -E '2~2s/^(.{8})(.{8}).{4}(.{6})/\1\2\3/g' > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi

    #remove adapter from Sci-Seq and swap barcode and UMI
    if [[ "$technology" == "sciseq" ]]; then
        for convFile in "${convFiles[@]}"; do
            #remove adapter if detected
            sed -E '
                /^ACGACGCTCTTCCGATCT/ {
                s/^(.{18})//g
                n
                n
                s/^(.{18})//g
                }'  $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
            #swap UMI and barcode
            sed -E '2~2s/(.{8})(.{10})/\2\1/' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile            
        done
    fi

    #remove adapter from SureCell (and correct phase blocks)
    if [[ "$technology" == "surecell" ]];
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC/ {
                s/.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC.*/\1\2\3\4/g
                n
                n
                s/.*(.{6}).{15}(.{6}).{15}(.{6}).{3}(.{8}).{3}/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi

    #converting barcodes
    echo " adjusting barcodes of R1 files"
    if [[ $barcodeadjust != 0 ]]; then
        if [[ $barcodeadjust -gt 0 ]]; then
            for convFile in "${convFiles[@]}"; do
                echo " handling $convFile ..."
                sed -i "2~2s/^.{${barcodeadjust}}//" $convFile #Trim the first n characters from the beginning of the sequence and quality
                echo "  ${convFile} adjusted"
           done
        elif [[ 0 -gt $barcodeadjust ]]; then
            for convFile in "${convFiles[@]}"; do
                echo " handling $convFile ..."
                toS=`printf '%0.sA' $(seq 1 $(($barcodeadjust * -1)))`
                toQ=`printf '%0.sI' $(seq 1 $(($barcodeadjust * -1)))`
                sed -i "2~4s/^/$toS/" $convFile #Trim the first n characters from the beginning of the sequence
                sed -i "4~4s/^/$toQ/" $convFile #Trim the first n characters from the beginning of the quality
                echo "  ${convFile} adjusted"
            done
        fi
    fi
    
    #UMI
    echo " adjusting UMIs of R1 files"
    if [[ 0 -gt $umiadjust ]]; then 
        for convFile in "${convFiles[@]}"; do
            echo " handling $convFile ..."
            toS=`printf '%0.sA' $(seq 1 $(($umiadjust * -1)))`
            toQ=`printf '%0.sI' $(seq 1 $(($umiadjust * -1)))`
            keeplength=`echo $((${barcode_default}+${umi_default}-($umiadjust * -1)))`
            sed -i "2~2s/^\(.\{${keeplength}\}\).*/\1/" $convFile #Trim off everything beyond what is needed
            sed -i "2~4s/$/$toS/" $convFile #Add n characters to the end of the sequence
            sed -i "4~4s/$/$toQ/" $convFile #Add n characters to the end of the quality
            echo "  ${convFile} adjusted"
        done
    fi
fi
##########



#####setting parameters for cellranger#####
d=""
if [[ -n $description ]]; then
    d="--description=$description"
fi
n=""
if [[ -n $ncells ]]; then
    n="--force-cells=$ncells"
fi
j=""
l=""
m=""
if [[ -n $jobmode ]]; then
    j="--jobmode=$jobmode"
    if [[ $jobmode == "local" ]]; then
        if [[ -n $ncores ]]; then
            l="--localcores=$ncores"
        fi
        if [[ -n $mem ]]; then
           m="--localmem=$mem"
        fi
    else
         if [[ -n $mem ]]; then
             m="--mempercore=$mem"
         fi
    fi
else
    if [[ -n $ncores ]]; then
        l="--localcores=$ncores"
    fi
    if [[ -n $mem ]]; then
        m="--localmem=$mem"
     fi
fi

#outputting command
echo "running cellranger ..."
echo ""
echo "#####cellranger command#####"

start=`date +%s`
echo "cellranger count --id=$id \
        --fastqs=$crIN \
        --lanes=$LANE \
        --r1-length=$totallength \
        --chemistry=$chemistry \
        --transcriptome=$reference \
        --sample=$SAMPLE \
        $d \
        $n \
        $j \
        $l \
        $m
"
echo "##########"
##########



#####running cellranger#####
cellranger count --id=$id \
        --fastqs=$crIN \
        --lanes=$LANE \
        --r1-length=$totallength \
        --chemistry=$chemistry \
        --transcriptome=$reference \
        --sample=$SAMPLE \
        $d \
        $n \
        $j \
        $l \
        $m
#        --noexit
#        --nopreflight
end=`date +%s`
runtime=$((end-start))
echo "cellranger run complete"
##########



#####remove files if convert is not running elsewhere#####
echo "updating .lock file"

#remove current job from counter (successfully completed)
lock=`cat ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.lock`
lock=$(($lock-1))
echo $lock > $lockfile

#check if jobs running
if [[ $lock -ge 1 ]]; then
    echo " total of $lock jobs for $lastcall technology are still run by cellranger ${cellrangerversion} in ${cellrangerpath}"
else
    echo " no other jobs currently run by cellranger ${cellrangerversion} in ${cellrangerpath}"
    echo " no conflicts: whitelist can now be changed for other technologies"
    rm -f $lockfile
fi
##########



#####Readjusting the barcodes in the cellranger output back to its original state##### 
if [[ ${barcodefile} != "default" ]]; then
    echo "replacing modified barcodes with the original in the output gene barcode matrix"
    perl ${BARCODERECOVER} ${barcodefile} ${barcodeadjust} ${id}
    echo "barcodes recovered"
fi
##########



#####extracting per cell data#####
if [[ $percelldata == true ]]; then
    echo "generating basic run statistics and per cell data"
    perl ${PERCELLSTATS} ${barcodefile} ${barcodeadjust} ${id}
    echo "per cell data generated"
fi
##########



#####printing out log#####
log="
#####Conversion tool log#####
cellranger ${cellrangerversion}

Original barcode format: ${technology} (then converted to 10x)

cellranger runtime: ${runtime}s
##########
"
echo "$log"
##########

exit 0

