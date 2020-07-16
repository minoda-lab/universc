#!/bin/sh

install=false

######convert version#####
convertversion="1.0.1"
##########



#####locate cellranger and get cellranger version#####
cellrangerpath=`which cellranger` #location of cellranger
if [[ -z $cellrangerpath ]]; then
    echo "cellranger command is not found."
    exit 1
fi
cellrangerversion=$(cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')')
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

TOOLS=${SDIR}/sub
BARCODERECOVER=${TOOLS}/RecoverBarcodes.pl
MAKEINDROPBARCODES=${TOOLS}/MakeIndropBarcodes.pl
PERCELLSTATS=${TOOLS}/ExtractBasicStats.pl
##########



#####define set options#####
lockfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.lock #path for .lock file
lastcallfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.last_called #path for .last_called
lastcall=`[[ -e $lastcallfile ]] &&  cat $lastcallfile || echo ""`
lastcall_b=`echo ${lastcall} | cut -f1 -d' '`
lastcall_u=`echo ${lastcall} | cut -f2 -d' '`
lastcall_p=`echo ${lastcall} | cut -f3 -d' '`
barcodedir=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes #folder within cellranger with the whitelist barcodes
barcodefile=""
crIN=input4cellranger #name of the directory with all FASTQ and index files given to cellranger
whitelistdir=${SDIR}/whitelists #path to whitelists
whitelistfile="outs/whitelist.txt" #name of the whitelist file added to the cellranger output
percellfile="outs/basic_stats.txt" #name of the file with the basic statistics of the run added to the cellranger output
##########



#####checki if convert and cellranger are writable#####
#cellranger
if ! [[ -w "$barcodedir" ]]; then
    echo "Error: Trying to run cellranger installed at ${cellrangerpath}"
    echo "launch_universc.sh can only run with cellranger installed locally"
    echo "Install cellranger in a directory with write permissions such as /home/`whoami`/local and export to the PATH"
    echo "The following versions of cellranger are found:"
    echo " `whereis cellranger`"
    exit 1
fi

#convert
if ! [[ -w "$SDIR" ]]; then
    echo "Error: Trying to run launch_universc.sh installed at $SDIR"
    echo "$SDIR must be writable to run launch_universc.sh"
    echo "Install launch_universc.sh in a directory with write permissions such as /home/`whoami`/local and export to the PATH"
    exit 1
fi
##########



#####usage statement#####
## detect shell across different OS
if [[ -z $VENDOR ]]; then
   if [[ $(uname -a | grep "[Mm]ac") ]]; then
       VENDOR="apple"
   else
       VENDOR="linux"
   fi
fi
if [[ $VENDOR != "apple" ]]; then
    SHELL=$(readlink -f /proc/$$/exe | cut -d'/' -f3)
else
    SHELL=$(ps -p $$ | awk '$1 == PP {print $4}' PP=$$)
fi
## detect how called (e.g., bash converrt.sh or ./launch_universc.sh)
if [[ $(which launch_universc.sh) == *"/"* ]] || [[ $0 == *"/"* ]]; then
    SHELL=''
    invocation=$0
else
    if [[ ! -z $ZSH_VERSION ]]; then
        SHELL="zsh"
    elif [[ ! -z $KSH_VERSION ]]; then
       SHELL="ksh"
    elif [[ ! -z $FISH_VERSION ]]; then
       SHELL="fish"
    elif [[ ! -z $BASH_VERSION ]]; then
        SHELL="bash"
    else
       SHELL=$SHELL
    fi
    invocation=$(echo $(basename $0))
fi
help='
Usage:
  '$SHELL' '$invocation' --testrun -t TECHNOLOGY
  '$SHELL' '$invocation' -t TECHNOLOGY --setup
  '$SHELL' '$invocation' -R1 FILE1 -R2 FILE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  '$SHELL' '$invocation' -R1 READ1_LANE1 READ1_LANE2 -R2 READ2_LANE1 READ2_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  '$SHELL' '$invocation' -f SAMPLE_LANE -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  '$SHELL' '$invocation' -f SAMPLE_LANE1 SAMPLE_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  '$SHELL' '$invocation' -v
  '$SHELL' '$invocation' -h

Convert sequencing data (FASTQ) from Nadia or ICELL8 platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
       --testrun                Initiates a test trun with the test dataset
  -R1, --read1 FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2 FILE             Read 2 FASTQ file to pass to cellranger
  -I1, --index1 FILE            Index (I1) FASTQ file to pass to cellranger (OPTIONAL)
  -I2, --index2 FILE            Index (I2) FASTQ file to pass to cellranger (OPTIONAL and EXPERIMENTAL)
  -f,  --file NAME              Path and the name of FASTQ files to pass to cellranger (prefix before R1 or R2)
                                  e.g. /path/to/files/Example_S1_L001
  
  -i,  --id ID                  A unique run id, used to name output folder
  -d,  --description TEXT       Sample description to embed in output files.
  -r,  --reference DIR          Path of directory containing 10x-compatible reference.
  -t,  --technology PLATFORM    Name of technology used to generate data.
                                Supported technologies:
                                  10x Genomics (version automatically detected): 10x, chromium
                                  10x Genomics version 2 (16bp barcode, 10bp UMI): 10x-v2, chromium-v2
                                  10x Genomics version 3 (16bp barcode, 12bp UMI): 10x-v3, chromium-v3
                                  CEL-Seq (8bp barcode, 4bp UMI): celseq
                                  CEL-Seq2 (6bp UMI, 6bp barcode): celseq2
                                  Drop-Seq (12bp barcode, 8bp UMI): nadia, dropseq
                                  ICELL8 version 3 (11bp barcode, 14bp UMI): icell8 or custom
                                  inDrops version 1 (19bp barcode, 6bp UMI): indrops-v1, 1cellbio-v1
                                  inDrops version 2 (19bp barcode, 6bp UMI): indrops-v2, 1cellbio-v2
                                  inDrops version 3 (8bp barcode, 6bp UMI): indrops-v3, 1cellbio-v3
                                  MARS-Seq (6bp barcode, 10bp UMI): marsseq, marsseq-v1
                                  MARS-Seq2 (7bp barcode, 8bp UMI): marsseq2, marsseq-v2
                                  Quartz-Seq2 (14bp barcode, 8bp UMI): quartzseq2-384
                                  Quartz-Seq2 (15bp barcode, 8bp UMI): quartzseq2-1536
                                  Sci-Seq (8bp UMI, 10bp barcode): sciseq
                                  SCRB-Seq (6bp barcode, 10bp UMI): scrbseq, mcscrbseq
                                  SeqWell (12bp barcode, 8bp UMI): seqwell
                                  Smart-seq2-UMI, Smart-seq3 (11bp barcode, 8bp UMI): smartseq
                                  SPLiT-Seq (10bp UMI, 18bp barcode): splitseq
                                  SureCell (18bp barcode, 8bp UMI): surecell, ddseq, biorad
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
#variable options
setup=false
convert=true
testrun=false
read1=()
read2=()
index1=()
index2=()
SAMPLE=""
LANE=()
id=""
description=""
reference=""
ncells=""
chemistry=""
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
        -I1|--index1|--index)
            shift
            if [[ "$1" != "" ]]; then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]; do
                    index1+=("${1/%\//}")
                    shift
                    arg=$1
                done
                next=true
            elif [[ -z $index1 ]]; then
                echo "Error: file input missing for --index1"
                exit 1
            fi
            ;;
        -I2|--index2)
            shift
            if [[ "$1" != "" ]]; then
                arg=$1
                while [[ ! "$arg" == "-"* ]] && [[ "$arg" != "" ]]; do
                    index2+=("${1/%\//}")
                    shift
                    arg=$1
                done
                next=true
            elif [[ -z $index1 ]]; then
                echo "Error: file input missing for --index1"
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



#####check if this is a test run#####
if [[ $testrun == "true" ]]; then
    if [[ ${#read1[@]} -gt 0 ]] || [[ ${#read2[@]} -gt 0 ]]; then
        echo "Error: for test run, no R1 or R2 file can be selected."
        exit 1
    elif [[ -n $reference ]]; then
        echo "Error: for test run, reference will be set automatically."
        exit 1
    elif [[ -n $id ]]; then
        echo "Error: for test run, id will be set automatically."
        exit 1
    fi
     
    TESTDIR=${SDIR}/test
    reference=${TESTDIR}/cellranger_reference/cellranger-tiny-ref/3.0.0
    id=test-tiny-${technology}
    description=${id}
    percelldata=true
    
    if [[ $technology == "10x" ]]; then
        gunzip -fk ${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L00[12]_R[12]_001.fastq.gz
        read1=("${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R1_001.fastq" "${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R1_001.fastq")
        read2=("${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_R2_001.fastq" "${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_R2_001.fastq")
        index1=("${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L001_I1_001.fastq.gz" "${TESTDIR}/shared/cellranger-tiny-fastq/3.0.0/tinygex_S1_L002_I1_001.fastq.gz")
    elif [[ $technology == "nadia" ]]; then
        gunzip -fk ${TESTDIR}/shared/dropseq-test/SRR1873277_S1_L001_R[12]_001.fastq.gz
        read1=("${TESTDIR}/shared/dropseq-test/SRR1873277_S1_L001_R1_001.fastq")
        read2=("${TESTDIR}/shared/dropseq-test/SRR1873277_S1_L001_R2_001.fastq")
    elif [[ $technology == "icell8" ]]; then
        gunzip -fk ${TESTDIR}/shared/icell8-test/ICELL8_01_S1_L00[12]_R[12]_001.fastq.gz
        read1=("${TESTDIR}/shared/icell8-test/ICELL8_01_S1_L001_R1_001.fastq" "${TESTDIR}/shared/icell8-test/ICELL8_01_S1_L002_R1_001.fastq")
        read2=("${TESTDIR}/shared/icell8-test/ICELL8_01_S1_L001_R2_001.fastq" "${TESTDIR}/shared/icell8-test/ICELL8_01_S1_L002_R2_001.fastq")
	barcodefile=${TESTDIR}/shared/icell8-test/BarcodeList.txt
	ncells=`wc -l $barcodefile | cut -f1 -d' '`
    else
        echo "Error: for test run, -t needs to be 10x, nadia, or icell8"
        exit 1
    fi
fi
##########



#####Check if technology matches expected inputs#####
if [[ $verbose ]]; then
    echo " checking option: --technology"
fi
if [[ "$technology" == "10x" ]] || [[ "$technology" == "chromium" ]]; then
    technology="10x"
elif [[ "$technology" == "10x-v2" ]] || [[ "$technology" == "chromium-v2" ]]; then
    technology="10x-v2"
elif [[ "$technology" == "10x-v3" ]] || [[ "$technology" == "chromium-v3" ]]; then
    technology="10x-v3"
elif [[ "$technology" == "celseq" ]] || [[ "$technology" == "cel-seq" ]]; then
    technology="celseq"
elif [[ "$technology" == "celseq2" ]] || [[ "$technology" == "cel-seq2" ]]; then
    technology="celseq2"
elif [[ "$technology" == "nadia" ]] || [[ "$technology" == "dropseq" ]] || [[ "$technology" == "drop-seq" ]]; then
    technology="nadia"
elif [[ "$technology" == "icell8" ]] || [[ "$technology" == "icell-8" ]]; then
    technology="icell8"
elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrops-v1" ]] || [[ "$technology" == "indropv1" ]] || [[ "$technology" == "indropsv1" ]] || [[ "$technology" == "1cellbio-v1" ]] || [[ "$technology" == "1cellbiov1" ]]; then
    technology="indrop-v1"
elif [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrops-v2" ]] || [[ "$technology" == "indropv2" ]] || [[ "$technology" == "indropsv2" ]] || [[ "$technology" == "1cellbio-v2" ]] || [[ "$technology" == "1cellbiov2" ]]; then
    technology="indrop-v2"    
elif [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "indrops-v3" ]] || [[ "$technology" == "indropv3" ]] || [[ "$technology" == "indropsv3" ]] || [[ "$technology" == "1cellbio-v3" ]] || [[ "$technology" == "1cellbiov3" ]]; then
    technology="indrop-v3"
elif [[ "$technology" == "marsseq" ]] || [[ "$technology" == "mars-seq" ]] || [[ "$technology" == "marsseq-v1" ]] || [[ "$technology" == "mars-seq-v1" ]] || [[ "$technology" == "marsseqv1" ]] || [[ "$technology" == "mars-seqv1" ]]; then
    technology="marsseq-v1"
elif [[ "$technology" == "marsseq2" ]] || [[ "$technology" == "mars-seq2" ]] || [[ "$technology" == "marsseq-v2" ]] || [[ "$technology" == "mars-seq-v2" ]] || [[ "$technology" == "marsseqv2" ]] || [[ "$technology" == "mars-seqv2" ]]; then
    technology="marsseq-v2"
elif [[ "$technology" == "quartz-seq2-384" ]] || [[ "$technology" == "quartzseq2-384" ]] || [[ "$technology" == "quartz-seq2-v3.1" ]] || [[ "$technology" == "quartzseq2-v3.1" ]] || [[ "$technology" == "quartzseq2v3.1" ]]; then
    technology="quartz-seq2-384"
elif [[ "$technology" == "quartz-seq2-1536" ]] || [[ "$technology" == "quartzseq2-1536" ]] || [[ "$technology" == "quartz-seq2-v3.2" ]] || [[ "$technology" == "quartzseq2-v3.2" ]] || [[ "$technology" == "quartzseq2v3.2" ]]; then
    technology="quartz-seq2-1536"
elif [[ "$technology" == "sciseq" ]] || [[ "$technology" == "sci-seq" ]]; then
    technology="sciseq"
elif [[ "$technology" == "scrbseq" ]] || [[ "$technology" == "scrb-seq" ]] || [[ "$technology" == "mcscrbseq" ]] || [[ "$technology" == "mcscrb-seq" ]]; then
    technology="scrbseq"
elif [[ "$technology" == "seqwell" ]] || [[ "$technology" == "seq-well" ]]; then
    technology="seqwell"
elif [[ "$technology" == "smartseq" ]] || [[ "$technology" == "smart-seq" ]] || [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smart-seq2" ]] ||  [[ "$technology" == "smartseq2-umi" ]] || [[ "$technology" == "smart-seq2-umi" ]] ||  [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "smart-seq3" ]]; then
    technology="smartseq"
elif [[ "$technology" == "splitseq" ]] || [[ "$technology" == "split-seq" ]]; then
    technology="splitseq"
elif [[ "$technology" == "surecell" ]] || [[ "$technology" == "surecellseq" ]] || [[ "$technology" == "surecell-seq" ]] || [[ "$technology" == "ddseq" ]] || [[ "$technology" == "dd-seq" ]] || [[ "$technology" == "bioraad" ]]; then
    technology="surecell"
elif [[ "$technology" == "custom"* ]]; then
    fields=$((`echo $technology | grep -o "_" | wc -l`+1))
    if [[ $fields -ne 3 ]]; then
        echo "Error: custom input must have exactly 3 fields separated by '_', e.g. custom_10_16"
    else
        b=`echo $technology | cut -f2 -d'_'`
        u=`echo $technology | cut -f3 -d'_'`
        if ! [[ "$b" =~ ^[0-9]+$ ]] || ! [[ "$u" =~ ^[0-9]+$ ]]; then
            echo "Error: option -t needs to be a technology listed or custom_<barcode>_<UMI>"
            exit 1
        fi
    fi
else
    echo "Error: option -t needs to be a technology listed or custom_<barcode>_<UMI>"
    exit 1
fi

if [[ $verbose ]]; then
    echo " technology set to ${technology}"
fi

if [[ "$technology" == "icell8" ]]; then
    echo "***WARNING: ${technology} should only be used for kits that have valid UMIs***"
fi
if [[ "$technology" == "smartseq" ]]; then
    echo "***WARNING: ${technology} should only be used for kits that have UMIs***"
fi
if [[ "$technology" == "smartseq" ]] || [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
    echo "***WARNING: launch_universc.sh does not support barcodes in dual indexes. Make sure that the R1 file is adjusted accordingly prior to running launch_universc.sh***"
fi
##########



#####Setting barcode and UMI length according to the set technology#####
if [[ $verbose ]]; then
    echo " setting barcode and UMI lengths for ${technology}"
fi

#barcode and umi lengths given by options
barcodelength=""
if [[ $verbose ]]; then
    echo "barcodelength before $barcodelength"
fi
umilength=""
minlength=""
if [[ "$technology" == "10x" ]]; then
    barcodelength=16
    umilength=10
    minlength=16
elif [[ "$technology" == "10x-v2" ]]; then
    barcodelength=16
    umilength=10
    minlength=16
elif [[ "$technology" == "10x-v3" ]]; then
    barcodelength=16
    umilength=12
    minlength=16
elif [[ "$technology" == "celseq" ]]; then
    barcodelength=8
    umilength=4
    minlength=8
elif [[ "$technology" == "celseq2" ]]; then
    barcodelength=6
    umilength=6
    minlength=6
elif [[ "$technology" == "nadia" ]]; then
    barcodelength=12
    umilength=8
    minlength=12
elif [[ "$technology" == "icell8" ]]; then
    barcodelength=11
    umilength=14
    minlength=11
elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
    barcodelength=19 
    umilength=6
    minlength=16
elif [[ "$technology" == "indrop-v3" ]]; then
    barcodelength=11
    umilength=6
    minlength=8
elif [[ "$technology" == "marsseq-v1" ]]; then
    barcodelength=6
    umilength=10
    minlength=6
elif [[ "$technology" == "marsseq-v2" ]]; then
    barcodelength=7
    umilength=8
    minlength=7
elif [[ "$technology" == "quartz-seq2-384" ]]; then
    barcodelength=14
    umilength=8
    minlength=14
elif [[ "$technology" == "quartz-seq2-1536" ]]; then
    barcodelength=15
    umilength=8
    minlength=15
elif [[ "$technology" == "sciseq" ]]; then
    barcodelength=10
    umilength=8
    minlength=10
elif [[ "$technology" == "scrbseq" ]]; then
    barcodelength=6 
    umilength=10
    minlength=6
elif [[ "$technology" == "seqwell" ]]; then
    barcodelength=8
    umilength=12
    minlength=8
elif [[ "$technology" == "smartseq" ]]; then
    barcodelength=11
    umilength=8
    minlength=11
elif [[ "$technology" == "splitseq" ]]; then
    barcodelength=18
    umilength=10
    minlength=18
elif [[ "$technology" == "surecell" ]]; then
    barcodelength=18
    umilength=8
    minlength=18
elif [[ "$technology" == "custom"* ]]; then
    barcodelength=`echo $technology | cut -f2 -d'_'`
    umilength=`echo $technology | cut -f3 -d'_'`
    minlength=${barcodelength}
fi

if [[ $minlength -gt 16 ]]; then
    minlength=16
fi

if [[ $verbose ]]; then
    echo " barcode and UMI lengths set as ${barcodelength} and ${umilength} respectively"
fi
##########



#####Setting chemisty#####
#check if chemistry matches expected inputs
if [[ $verbose ]]; then
    echo " checking option: --chemistry"
fi

if [[ -n $chemistry ]]; then
    if [[ "$chemistry" != "auto" ]] && \
    [[ "$chemistry" != "threeprime" ]] && \
    [[ "$chemistry" != "fiveprime" ]] && \
    [[ "$chemistry" != "SC3Pv1" ]] && \
    [[ "$chemistry" != "SC3Pv2" ]] && \
    [[ "$chemistry" != "SC3Pv3" ]] && \
    [[ "$chemistry" != "SC5P-PE" ]] && \
    [[ "$chemistry" != "SC5P-R2" ]] && \
    [[ "$chemistry" != "SC-FB" ]]; then
        echo "Error: option -c can be auto, threeprime, fiveprime, SC3Pv1, SC3Pv2, SC3Pv3, SC5P-PE, SC5P-R2, or SC-FB"
	exit 1
    fi
fi

#determine what chemistry is recommended
temp_chemistry="SC3Pv2"
if [[ $umilength -gt 10 ]]; then
    temp_chemistry="SC3Pv3"
fi
if [[ "$technology" == "10x" ]]; then
    temp_chemistry="auto"
fi

#set chemistry
if [[ -z ${chemistry} ]]; then
    chemistry=${temp_chemistry}
elif [[ "$chemistry" != "$temp_chemistry" ]]; then
    echo "***WARNING: chemistry is set to ${chemistry} where ${temp_chemistry} would have been chosen automatically. proceed with causion.***"
fi

#set default barcode and umi lengths
barcode_default=16
umi_default=10
if [[ "$chemistry" == "SC3Pv3" ]]; then
    umi_default=12
fi
totallength=`echo $((${barcode_default}+${umi_default}))`

#adjustment lengths
if [[ $verbose ]]; then
   echo "barcode length: $barcodelength"
   echo "barcode default: $barcode_default"
fi
barcodeadjust=`echo $(($barcodelength-$barcode_default))`
umiadjust=`echo $(($umilength-$umi_default))`
if [[ $verbose ]]; then
   echo "barcode adjust: $barcodeadjust "
fi

if [[ $verbose ]]; then
    echo " chemistry set to ${chemistry}"
fi
##########



#####Input file curation 1: Check if number of R1, R2, I1, and I2 files are as expected#####
#get unique inputs
read1=(`echo "${read1[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '`)
read2=(`echo "${read2[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '`)
index1=(`echo "${index1[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '`)
index2=(`echo "${index2[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '`)

#check for presence of read1 and read2 files
if [[ $verbose ]]; then
    echo " checking option: --read1 and --read2"
fi

if [[ $setup == "false" ]]; then
    if [[ ${#read1[@]} -eq 0 ]]; then
        echo "Error: option --read1 or --file is required"
        exit 1
    elif [[ ${#read2[@]} -eq 0 ]]; then
        echo "Error: option --read2 or --file is required"
        exit 1
    fi
fi

#check for presence of indexe files
if [[ $verbose ]]; then
    echo " checking option: --index1 and --index2"
fi

#index 1
if [[ $setup == "false" ]]; then
    if [[ ${#index1[@]} -ne ${#read1[@]} ]]; then
        if [[ ${#index1[@]} -gt 0 ]]; then
            echo " Error: number of index1 files is not matching the number of R1 and R2 files"
            echo " either give no index1 file, or give index1 file for each and every read1 file"
            exit 1
        fi
    fi
fi

#index 2
if [[ $setup == "false" ]]; then
    #only check I2 for dual-indexed techniques
    if [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "sci-seq" ]] || [[ "$technology" == "smartseq" ]]; then
        if [[ ${#index2[@]} -ne ${#read1[@]} ]]; then        
            echo " Error: number of index1 files is not matching the number of index2 files"
            echo " for $technology, either give no index files or give index1 and index2 for each and every read1 and read2 file"
            exit 1
        fi
    elif [[ ${#index2[@]} -gt 0 ]]; then
        echo " Error: $technology does not support dual index"
        echo " re-run without selecting any files for --index2"
        exit 1
    fi   
fi

if [[ $verbose ]]; then
    echo "  ${#read1[@]} read1s: ${read1[@]}"
    echo "  ${#read2[@]} read2s: ${read2[@]}"
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]} I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]} I2s: ${index2[@]}"
    fi
    echo "  number of these files are as expected"
fi

keys=("R1" "R2")
if [[ ${#index2[@]} -gt 0 ]]; then
    keys=("R1" "R2" "I1" "I2")
elif [[ ${#index1[@]} -gt 0 ]]; then
    keys=("R1" "R2" "I1")
fi
##########



#####Input file curation 2: Check R1, R2, I1, and I2 files for their extensions#####
##allows incomplete file names and processing compressed files
if [[ $verbose ]]; then
    echo "key: ${keys[@]}"
fi
for key in ${keys[@]}; do
    readkey=$key
    list=()
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    elif [[ $readkey == "I1" ]]; then
        list=("${index1[@]}")
    elif [[ $readkey == "I2" ]]; then
        list=("${index2[@]}")
    fi
    
    for j in ${!list[@]}; do
        read=${list[$j]}
        if [[ $verbose ]]; then
            echo "  checking file format for $read ..."
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
        
        if [[ $verbose ]]; then
             echo "  $read processed"
        fi
        
        list[$j]=$read
    done
    
    if [[ $readkey == "R1" ]]; then
        read1=("${list[@]}")
    elif [[ $readkey == "R2" ]]; then
        read2=("${list[@]}")
    elif [[ $readkey == "I1" ]]; then
        index1=("${list[@]}")
     elif [[ $readkey == "I2" ]]; then
        index2=("${list[@]}")
    fi
done

if [[ $verbose ]]; then
    echo "  ${#read1[@]} read1s: ${read1[@]}"
    echo "  ${#read2[@]} read2s: ${read2[@]}"
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]} I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]} I2s: ${index2[@]}"
    fi
    echo "  files exist, with extentions compatible with launch_universc.sh"
fi
##########



#####Input file curation 3: renaming read1, read2, index1, and index2 file name if not compatible with the launch_universc.sh#####
for key in ${keys[@]}; do
    readkey=$keys
    list=""
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    elif [[ $readkey == "I1" ]]; then
         list=("${index1[@]}")
    elif [[ $readkey == "I2" ]]; then
         list=("${index2[@]}")
    fi
    
    for j in ${!list[@]}; do
        read=${list[$j]}
        if [[ $verbose ]]; then
            echo "  checking file name for $read ..."
        fi
        
        if [[ -h $read ]]; then
            path=`readlink -f $read`
            if [[ $verbose ]]; then
                echo "***Warning: file $read not in current directory. Path to the file captured instead.***"
                echo " (file) $read"
                echo " (path) $path"
            fi
            read=${path}
        fi
        case $read in
            #check if contains lane before read
            *_L0[0123456789][0123456789]_$readkey*)
                if [[ $verbose ]]; then
                    echo "  $read compatible with lane"
                fi
            ;;
            *) 
                #rename file
                if [[ $verbose ]]; then
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
                if [[ $verbose ]]; then
                    echo "  $read compatible with sample"
                fi
            ;;
            *)
                #rename file
                if [[ $verbose ]]; then
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
                if [[ $verbose ]]; then
                    echo "  $read compatible with suffix"
                fi
            ;;
            *)
                #rename file
                if [[ $verbose ]]; then
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
     elif [[ $readkey == "I1" ]]; then
        index1=("${list[@]}")
      elif [[ $readkey == "I2" ]]; then
        index2=("${list[@]}")
    fi
done

if [[ $verbose ]]; then
    echo "  ${#read1[@]} read1s: ${read1[@]}"
    echo "  ${#read2[@]} read2s: ${read2[@]}"
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]} I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]} I2s: ${index2[@]}"
    fi
    echo "  names of these files are compatible with launch_universc.sh"
fi
##########



#####Input file curation 4: Technology-specific adjustments#####
#reorganize for indrop-v3
if [[ "$technology" == "indrop-v3" ]]; then
    #check for indexes in R2 and R3 (R1 -> R2; R2 -> I1; R3 -> I2; R4 -> R1)
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo " index I1 found for ${technology}"
    else
        #checking for R2 and R3 index files
        echo " checking for index I1 in R2 files..."
        for ii in $(seq 1 1 ${#read1[@]}); do
             #iterate over read1 inputs
             indexfile=${read1[$(( $ii -1 ))]}
             #derive I1 filename from R1 filename
             indexfile=$(echo $indexfile | perl -pne 's/(.*)_R1/$1_R2/' )
             #only add index2 files to list variable if file exists (this will only run if I1 not found above)
             if [[ -f $indexfile ]] || [[ -f ${indexfile}.gz ]] || [[ -f $indexfile.fastq ]] || [[ -f ${indexfile}.fastq.gz ]] || [[ -f $indexfile.fq ]] || [[ -f ${indexfile}.fq.gz ]]; then
                 index1+=("$indexfile")
             fi
        done
    fi
    if [[ ${#index2[@]} -eq ${#read1[@]} ]] && [[ ${#index2[@]} -ge 1 ]]; then
        echo " index I2 found for ${technology}"
    else
        #checking for R2 and R3 index files
        echo " checking for index I2 in R3 files..."
        for ii in $(seq 1 1 ${#read1[@]}); do  
             #iterate over read1 inputs
             indexfile=${read1[$(( $ii -1 ))]}
             #derive I2 filename from R1 filename
             indexfile=$(echo $indexfile | perl -pne 's/(.*)_R1/$1_R3/' )
             #only add index2 files to list variable if file exists (this will only run if I2 not found above)
             if [[ -f $indexfile ]] || [[ -f ${indexfile}.gz ]] || [[ -f $indexfile.fastq ]] || [[ -f ${indexfile}.fastq.gz ]] || [[ -f $indexfile.fq ]] || [[ -f ${indexfile}.fq.gz ]]; then
                 index2+=("$indexfile")
             fi
        done
    fi
    #checking for R1 and R4 index files
    echo " checking for read R4 files..."
    for ii in $(seq 1 1 ${#read1[@]}); do
         #iterate over read1 inputs
         indexfile=${read1[$(( $ii -1 ))]}
         #derive I2 filename from R1 filename
         indexfile=$(echo $indexfile | perl -pne 's/(.*)_R1/$1_R4/' )
         #only replace R2 files with R4 variable if R4 file exists (otherwise keep R2)
         if [[ -f $indexfile ]] || [[ -f ${indexfile}.gz ]] || [[ -f $indexfile.fastq ]] || [[ -f ${indexfile}.fastq.gz ]] || [[ -f $indexfile.fq ]] || [[ -f ${indexfile}.fq.gz ]]; then
             read2[$(( $ii -1 ))]=("$indexfile")
         fi
         #Note that R2 and R1 are inverted below (thus [I1, I2]R2 or [R2, R3]R4 are converted to R1)
    done
    #checking inDrops-v3 indexes are accepted (by checking that index 2 is correct)
    ##note that index 1 will be detected as I1 OR R2
    ##note that read1 will be detected R2 or R4
    if [[ ${#index2[@]} -eq ${#read1[@]} ]] && [[ ${#index2[@]} -ge 1 ]]; then
        echo " indexes ${index1[@]} and ${index2[@]} found for ${technology}"
    else
        if [[ $setup == "false" ]]; then
            echo "ERROR: note that ${technology} expects dual indexes: I1 and I2 OR R2 and R3"
            exit 1
        fi
    fi
fi

#inverting R1 and R2 for specific technologies
if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "splitseq" ]]; then
    #invert read1 and read2
    echo "***WARNING: technology is set to ${technology}. barcodes on Read 2 will be used***"
    tmp=$read1
    read1=$read2
    read2=$tmp
    tmp=""
fi

if [[ $verbose ]]; then
    echo "  ${#read1[@]} read1s: ${read1[@]}"
    echo "  ${#read2[@]} read2s: ${read2[@]}"
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]} I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]} I2s: ${index2[@]}"
    fi
    echo "  input files adjusted for technology-specific conditions"
fi
##########



#####Input file curation 5: Catpuring sample name#####
#checking the quality of fastq file names
read12=("${read1[@]}" "${read2[@]}")
if [[ ${#index2[@]} -gt 0 ]]; then
    read12=("${read1[@]}" "${read2[@]}" "${index1[@]}" "${index2[@]}")
elif [[ ${#index1[@]} -gt 0 ]]; then
    read12=("${read1[@]}" "${read2[@]}" "${index1[@]}")
fi

if [[ $verbose ]]; then
    echo "reads:" ${read12[@]}
fi

for fq in "${read12[@]}"; do
    if [[ $verbose ]]; then
        echo "read: $fq"
    fi
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

    if [[ $verbose == true ]]; then
       echo "SAMPLE: $SAMPLE"
       echo "sample (field): $sn"
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

if [[ $verbose ]]; then
     echo "  read1, read2, index1, and index2 file curation complete"
fi

LANE=$(echo "${LANE[@]}" | tr ' ' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
##########
 


#####Set the input barcode file#####
if [[ $verbose ]]; then
    echo " setting whitelist barcode file."
fi

if [[ -n "$barcodefile" ]]; then
    if [[ ! -f $barcodefile ]]; then
        echo "Error: File selected for option --barcodefile does not exist"
	exit 1
    else
        #getting absolute path
        barcodefile=$(readlink -f $barcodefile)
        #allowing WellList from ICELL8 and other well-based techniques
        if [[ "$technology" == "icell8" ]] || [[ "$technology" == "quartz-seq2*" ]] || [[ "$technology" == "splitseq" ]] || [[ "$technology" == "smartseq*" ]] || [[ "$technology" == "seqwell" ]] || [[ "$technology" == "sciseq" ]] || [[ "$technology" == "custom" ]]; then
            seg=$'\t'
            n_col=$(awk -F'\t' '{print NF}' $barcodefile | sort -nu | tail -n 1)
            if [[ $n_col -eq 1 ]]; then
                seg=","
                n_col=$(awk -F',' '{print NF}' $barcodefile | sort -nu | tail -n 1)
            fi
            
            if [[ $n_col -gt 1 ]]; then
                new_barcodefile=${barcodefile%.*}_barcode.txt
                
                #get column with barcodes
		col_n=$(head -n 1 $barcodefile | tr "${seg}" "\n" | grep -n "[Bb]arcode" | head -n 1 | cut -d":" -f1)
		if [[ -z $col_n ]]; then
                    col_n=1
                fi
                
                if [[ $(head -n 1 $barcodefile | grep "[Bb]arcode") ]]; then
                    #removes header (1st line) containing colname "barcode"
		    tail -n $(($(wc -l $barcodefile | cut -d" " -f1)-1)) $barcodefile | cut -f$col_n -d "$seg" > ${new_barcodefile}
                else
                    echo "***WARNING: barcode file has multiple columns and none named 'barcode(s)'"
                    echo "Note: please check that column 1 of the barcode file contains the barcodes as required"
                    head $barcodefile | cut -f1 -d "${seg}"
                    #assumes no headers
                    cut -f1 -d "${seg}" $barcodefile > ${new_barcodefile}
                fi
                barcodefile=${new_barcodefile}
            fi
        fi
    fi
else
    if [[ "$technology" == "10x"* ]]; then
        barcodefile="default:10x"
    elif [[ "$technology" == "icell8" ]]; then
        barcodefile=${whitelistdir}/ICELL8_barcode.txt
	echo "***WARNING: selected barcode file (${barcodefile}) contains barcodes for all wells in ICELL8. valid barcode will be an overestimate***"
    elif [[ "$technology" == "marsseq-v2" ]]; then
        barcodefile=${whitelistdir}/MARS-Seq2_barcode.txt
    elif [[ "$technology" == "quartz-seq2-384" ]]; then
        barcodefile=${whitelistdir}/Quartz-Seq2-384_barcode.txt
    elif [[ "$technology" == "quartz-seq2-1536" ]]; then
        barcodefile=${whitelistdir}/Quartz-Seq2-1536_barcode.txt
    elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        echo "***WARNING: whitelist for ${technology} is modified from the original barcodes (https://github.com/indrops/indrops/tree/master/ref/barcode_lists), first 8bp of list1 and list2 are joind to generate a 16bp barcode***"
        barcodelength=${minlength}
        if [[ $verbose ]]; then
            echo "  barcode adjusted to ${barcodelength}bp to match the length in the default whitelist for ${technology}"
        fi
        if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
            barcodefile=${whitelistdir}/inDrop-v2_barcodes.txt
        elif [[ "$technology" == "indrop-v3" ]]; then
            barcodefile=${whitelistdir}/inDrop-v3_barcodes.txt
            echo "***WARNING: ***combination of list1 and list2 from indrop-v2 (https://github.com/indrops/indrops/issues/32)***"  
        fi
    elif [[ "$technology" == "smartseq" ]]; then
        barcodefile=${whitelistdir}/SmartSeq3_barcode.txt 
    else
        echo "***WARNING: whitelist for ${technology} will be all possible combinations of ${minlength}bp. valid barcode will be 100% as a result***"
        barcodelength=${minlength}
	barcodefile=${whitelistdir}/AllPossibilities_${barcodelength}_barcodes.txt
    fi
fi

if [[ $verbose ]]; then
    echo " barcode file set to ${barcodefile}"
fi
##########



#####Generate whitelist file (if absent)#####
if [[ $verbose ]]; then
    echo " generating the selected barcode file"
fi
if [[ -f ${barcodefile} ]]; then
    if [[ $verbose ]]; then
        echo "  barcodefile file exists"
    fi
    if ! [[ "${barcodefile}" == "${whitelistdir}"* ]]; then
        if [[ $verbose ]]; then
            echo "  ensuring all barcode within selected whitelist are in upper case"
        fi
        sed -i 's/.*/\U&/g' $barcodefile
    fi
else
    if [[ ${barcodefile} == "default:10x" ]]; then
        if [[ $verbose ]]; then
            echo "  default 10x barcode whitelist will be used"
        fi
    elif ! [[ "${barcodefile}" == "${whitelistdir}"* ]]; then
        echo "Error: user selected barcode file (${barcodefile}) does not exist"
        exit 1
    else
        if [[ $verbose ]]; then
            echo "  generating a new barcode whitelist for ${technology}"
        fi
        if [[ "$technology" == "indrop-v"* ]]; then
            if [[ "$technology" == "indrop-v1" ]] || [[ $technology"" == "indrop-v2" ]]; then
                ${MAKEINDROPBARCODES} ${whitelistdir}/inDrop_gel_barcode1_list.txt ${whitelistdir}/inDrop_gel_barcode2_list.txt v2 ${whitelistdir}
            elif [[ "$technology" == "indrop-v3" ]]; then
                #ignore barcodes in index (R1 only)
                cp ${whitelistdir}/inDrop_gel_barcode2_list.txt ${whitelistdir}/inDrop-v3_barcodes.txt
                #${MAKEINDROPBARCODES} ${whitelistdir}/inDrop_gel_barcode1_list.txt ${whitelistdir}/inDrop_gel_barcode2_list.txt v3 ${whitelistdir}
            fi
        else
            #generating permutations of ATCG of barcode length (non-standard evaluation required to run in script)
            echo $(eval echo $(for ii in $(eval echo {1..${barcodelength}}); do echo "{A,T,C,G}"; done | tr "\n" " " | sed "s/ //g" | xargs -I {} echo {})) | sed 's/ /\n/g' | sort | uniq > ${barcodefile}
        fi
    fi
fi

if [[ $verbose ]]; then
    echo " barcodefile generated"
fi
##########



#####check if reference is present#####
if [[ -z $reference ]]; then
    if [[ $setup == "false" ]] || [[ ${#read1[@]} -ne 0 ]] || [[ ${#read2[@]} -ne 0 ]]; then
        echo "Error: option --reference is required";
        exit 1
    fi
fi
##########



#####check if ncells is an integer#####
int='^[0-9]+$'
if [[ -z "$ncells" ]]; then
    ncells=""
elif ! [[ $ncells =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --force-cells must be an integer"
    exit 1
fi
##########



#####check if ncores is an integer#####
int='^[0-9]+$'
if [[ -z "$ncores" ]]; then
    ncores=""
elif ! [[ $ncores =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --localcores must be an integer"
    exit 1
fi
##########



#####check if mem is a number#####
int='^[0-9]+([.][0-9]+)?$'
if [[ -z "$mem" ]]; then
    mem=""
elif ! [[ $mem =~ $int ]] && [[ $setup == "false" ]]; then
    echo "Error: option --localmem or --mempercore must be a number (of GB)"
    exit 1
fi
##########



#####check if chemistry matches expected input#####
#allow "auto" only for 10x
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
##########



#####checking if jobmode matches expected input#####
if [[ "$jobmode" != "local" ]] && [[ "$jobmode" != "sge" ]] && [[ "$jobmode" != "lsf" ]] && [[ "$jobmode" != *"template" ]]; then
    echo "Error: option --jobmode must be local, sge, lsf, or a .template file"
    exit 1
fi
##########



#####check if ID is present#####
if [[ -z $id ]]; then
    if [[ ${#read1[@]} -ne 0 ]] || [[ ${#read2[@]} -ne 0 ]]; then
        echo "Error: option --id is required"
        exit 1
    fi
fi
crIN=${crIN}_${id}
##########


#####checking if crIN exists#####
if [[ ! -d $crIN ]]; then
    convert=true
    echo "***Warning: conversion was turned on because directory $crIN was not found***"
fi
##########



#####check if UniverSC is running already#####
echo " checking if UniverSC is running already"

#set up .lock file
if [[ ! -f $lockfile ]]; then
    echo "  creating .lock file"
    echo 0 > $lockfile
    lock=`cat $lockfile`
else
    #check if jobs are running (check value in .lock file)
    echo "  checking .lock file"
    lock=`cat $lockfile`
    
    if [[ $lock -le 0 ]]; then
        echo "  call accepted: no other cellranger jobs running"
        lock=1
        if [[ $setup == "false" ]]; then 
                    echo $lock > $lockfile
        fi
    else
        if [[ -f $lastcallfile ]]; then
	    echo "  total of $lock cellranger ${cellrangerversion} jobs are already running in ${cellrangerpath} with barcode length (${lastcall_b}), UMI length (${lastcall_u}), and whitelist barcodes (${lastcall_p})"
            
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
echo "CONVERSION: $convert"
if [[ $convert == "false" ]]; then
    echo "***Warning: adjustment for barcode and UMI length was skipped***"
fi
echo "##########"
echo ""
##########



####setup whitelist#####
if [[ $verbose ]]; then
    echo "lock $lock"
fi
if [[ $lock -eq 0 ]]; then
    echo "setup begin"
    echo "updating barcodes in $barcodedir for cellranger version ${cellrangerversion} installed in ${cellrangerpath} ..."
    
    cd $barcodedir
    
    #restore assert functions if cellranger version is 3 or greater
    echo " restoring cellranger"
    if [[ $verbose  ]]; then
        echo "${cellrangerversion}"
        echo "${cellrangerversion} 3.0.0" | tr " " "\n" | sort -V | head -n 1
        echo "last call: $lastcall_p"
    fi
    if [[ $(echo "${cellrangerversion} 3.0.0" | tr " " "\n" | sort -V | tail -n 1)  == ${cellrangerversion} ]]; then
        if [[ $verbose  ]]; then
            echo "cellranger version 3.0.0 or greater"
        fi
        if [[ ! -z $barcodefile ]] && [[ $verbose  ]]; then
            echo "barcodefile: $barcodefile"
        fi
        if [[ $technology == "10x" ]] && [[ $barcodefile == "default:10x" ]]; then
            #restore checking barcodes
            if [[ $verbose  ]]; then
                echo "restore barcode checks"
            fi
            sed -i "s/#*#if gem_group == prev_gem_group/if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#*#assert barcode_idx >= prev_barcode_idx/assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#*#assert np.array_equal(in_mc.get_barcodes(), barcodes)/assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        elif [[ $lastcall_p == "default:10x" ]] || [[ ! -f $lastcallfile ]]; then
            #disable checking barcodes
            if [[ $verbose  ]]; then
                 echo "disable barcode checks"
            fi
            sed -i "s/if gem_group == prev_gem_group/#if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        fi
    fi
    if [[ $technology == "10x" ]] && [[ $barcodefile == "default:10x" ]]; then
        #restore cloupe generation
        if [[ $verbose  ]]; then
            echo "restore cloupe" 
        fi
        ## list cloupe output as (not null)
        if [[ $verbose  ]]; then
            echo "sed -i '/cloupe/s/null/CLOUPE_PREPROCESS\.output_for_cloupe/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i '/cloupe/s/null/CLOUPE_PREPROCESS\.output_for_cloupe/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro 
        ## add cloupe to outputs
        if [[ $verbose  ]]; then
            echo "sed -i '/out cloupe *cloupe/ {s/^#*#//g}' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i '/out cloupe *cloupe/ {s/^#*#//g}' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        #restore 11 lines for cloupe preprocess call (all following steps are needed to be called together or the call will break)
        ## restore defining CLOUPE_PREPROCESS
        if [[ $verbose  ]]; then
            echo "sed -i 's/^#*#@include "_cloupe_stages.mro"/@include "_cloupe_stages.mro"/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i 's/^#*#@include "_cloupe_stages.mro"/@include "_cloupe_stages.mro"/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        ## remove listing CLOUPE in output
        sed -i '/output_for_cloupe/s/^#*#//g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro 
        ## remove calling CLOUPE_PREPROCESS
        ### iterate over all files calling CLOUPE_PREPROCESS
        for file in $(grep -l  "call CLOUPE_PREPROCESS"  ${cellrangerpath}-cs/${cellrangerversion}/mro/*.mro )
        do
            #find  start of CLOUPE_PREPROCESS call
            num=$(grep -n "call CLOUPE_PREPROCESS" $file |  head -n 1 | cut -d":" -f1)
            #find end of CLOUPE_PREPROCESS call
            num2=$(($num +$(tail -n $(echo $(($(wc -l $file | cut -f1 -d " ") - $num))) $file | grep -n ")" | head -n 1 | cut -d":" -f1)))
            if [[ $verbose ]]; then
               echo "lines ${num}-${num2} restored in $file"
            fi
            eval "sed -i '$(echo "${num},${num2}s/^#*#//g")' $file"
        done
    elif [[ $lastcall_p == "default:10x" ]] || [[ ! -f $lastcallfile ]]; then
        #disable cloupe generation
        if [[ $verbose  ]]; then
            echo "disable cloupe"
        fi
        ## remove cloupe from outputs
        sed -i '/out cloupe *cloupe/ {s/^/#/g}' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        #remove 11 lines for cloupe preprocess call (all following steps are needed to be suppressed together or call will break)
        ## remove defining CLOUPE_PREPROCESS
        sed -i 's/@include "_cloupe_stages.mro"/#@include "_cloupe_stages.mro"/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        ## remove listing CLOUPE in output
        sed -i '/output_for_cloupe/s/^/#/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro 
         ## list cloupe output as null
        sed -i '/output_for_cloupe/s/CLOUPE_PREPROCESS\.output_for_cloupe/null/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        ## remove calling CLOUPE_PREPROCESS
        ### iterate over all files calling CLOUPE_PREPROCESS
        for file in $(grep -l  "call CLOUPE_PREPROCESS"  ${cellrangerpath}-cs/${cellrangerversion}/mro/*.mro ); do
            #find  start of CLOUPE_PREPROCESS call
            num=$(grep -n "call CLOUPE_PREPROCESS" $file |  head -n 1 | cut -d":" -f1)
            #find end of CLOUPE_PREPROCESS call
            num2=$(($num +$(tail -n $(echo $(($(wc -l $file | cut -f1 -d " ") - $num))) $file | grep -n ")" | head -n 1 | cut -d":" -f1)))
            if [[ $verbose ]]; then
               echo "lines ${num}-${num2} commented out of $file"
            fi
            eval "sed -i '$(echo "${num},${num2}s/^/#/g")' $file"
        done
    fi
    echo " ${cellrangerpath} set for $technology"
    
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
    
    #convert whitelist to the appropriate barcode
    echo " converting whitelist"
    if [[ ${barcodefile} == "default:10x" ]]; then
        #for version 2
        cp ${v2}.backup ${v2}
        #for version 3
        cp ${v3}.backup.gz ${v3}.gz
    else
        #for version 2
        cat ${barcodefile} > ${v2}
        echo "barcode adjust: $barcodeadjust "
        if [[ $barcodeadjust -gt 0 ]]; then
            sed -i "s/^.{${barcodeadjust}}//" ${v2} #Trim the first n characters from the beginning of the sequence and quality
        elif [[ 0 -gt $barcodeadjust ]]; then
            As=`printf '%0.sA' $(seq 1 $(($barcodeadjust * -1)))`
            echo As: $As
            sed -i "s/^/$As/" ${v2} #Trim the first n characters from the beginning of the quality
        fi
        #for version 3
        cat ${v2} > ${v3}
        gzip -f ${v3}
        rm translation/${v3}.gz
        ln -s ${v3}.gz translation/${v3}.gz
    fi
    echo " whitelist converted"

    echo "verbose $verbose"
    #change last call file
    if [[ $verbose ]]; then
        echo "Setting last call as..."
        echo "${barcodelength} ${umilength} ${barcodefile}"
    fi
    echo "${barcodelength} ${umilength} ${barcodefile}" > $lastcallfile

    cd - > /dev/null

    echo "setup complete"
fi
#########



#####exiting when setup is all that is requested#####
if [[ $setup == "true" ]]; then
    lock=`cat $lockfile`
    if [[ $lock -ge 1 ]]; then
        lock=$(($lock-1))
        echo $lock > $lockfile
    else
       rm -rf $lockfile
    fi
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

if [[ $verbose ]]; then
    echo "convert: $convert"
fi
if [[ $convert == "true" ]]; then
    echo "moving file to new location"
fi

if [[ $verbose ]]; then
    echo "Processing Read1"
    echo "Fastqs: ${read1[@]}"
    echo "${read1[@]}"
fi

crR1s=()
for fq in "${read1[@]}"; do
    if [[ $verbose ]]; then
        echo $fq
    fi
    to=`basename $fq`
    to="${crIN}/${to}"
    
    #invert read1 and read2
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        #where converted "read1" is R2 in source files (corrected names for cellranger)
        echo "using transcripts in Read 2 for ${technology}"
        to=`echo $to | sed -e "s/_R2_/_R1_/g"`
    fi
    
    if [[ $verbose ]]; then
        echo $to
    fi
    crR1s+=($to)
    
    echo " handling $fq ..."
    if [[ ! -f $to ]] || [[ $convert == "true" ]]; then
        if [[ $verbose ]]; then
            echo "cp -f $fq $to"
        fi
        cp -f $fq $to
    fi
    if [[ $convert == "true" ]]; then
        convFiles+=($to)
    fi
done

crR2s=()
if [[ $verbose ]]; then
     echo "Processing Read2"
     echo "Fastqs: ${read2[@]}"
     echo "${read2[@]}"
fi
for fq in "${read2[@]}"; do
    if [[ $verbose ]]; then
        echo "$fq";
    fi
    to=`basename $fq`
    to="${crIN}/${to}"
    to=$(echo "$to" | sed 's/\.gz$//')
    
    #invert read1 and read2
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        #where converted "read2" is R1 in source files
        echo "using transcripts in Read 1 for ${technology}"
        to=`echo $to | sed -e "s/_R1_/_R2_/g"`
        #where converted "read2" is R4 in source files
        if [[ "$technology" == "indrop-v3" ]]; then
            to=`echo $to | sed -e "s/_R4_/_R2_/g"`
        fi
    fi

    if [[ $verbose ]]; then
        echo "$to"
    fi
    crR2s+=($to)
    
    echo " handling $fq ..."
    if [[ ! -f $to ]] || [[ $convert == "true" ]]; then
        if [[ $verbose ]]; then 
            echo "cp -f $fq $to"; 
        fi
        cp -f $fq $to
    fi
done

if [[ ${#index1[@]} -ge 1 ]]; then
    crI1s=()
    if [[ $verbose ]]; then
         echo "Processing Index"
         echo "Fastqs: ${index1[@]}"
         echo "${index1[@]}"
    fi
    for fq in "${index1[@]}"; do
        if [[ $verbose ]]; then
            echo "$fq"
        fi
        to=`basename $fq`
        to="${crIN}/${to}"
        to=$(echo "$to" | sed 's/\.gz$//')
        
        #convert index for R2 and R3
        if [[ "$technology" == "indrop-v3" ]]; then
            #where converted "index1" is R2 in source files (corrected names for cellranger)
            echo "using transcripts in Read 2 for ${technology}"
            to=`echo $to | sed -e "s/_R2_/_I1_/g"`
        fi    
        
        if [[ $verbose ]]; then
            echo "$to"
        fi
        crI1s+=($to)
        
        echo "handling $fq ..."
        if [[ ! -f $to ]] || [[ $convert  == true ]]; then
            if [[ $verbose ]]; then
                echo "cp -f $fq $to"
	    fi
            cp -f $fq $to
        fi
    done
fi

if [[ ${#index2[@]} -ge 1 ]]; then
    crI2s=()
    if [[ $verbose ]]; then
        echo "Processing Index"
        echo "Fastqs: ${index2[@]}"
        echo "${index2[@]}"
    fi
    for fq in "${index2[@]}"; do
        if [[ $verbose ]]; then
            echo "$fq"
        fi
        to=`basename $fq`
        to="${crIN}/${to}"
        to=$(echo "$to" | sed 's/\.gz$//')
        
        #convert index for R2 and R3
        if [[ "$technology" == "indrop-v3" ]]; then
            #where converted "index2" is R3 in source files (corrected names for cellranger)
            echo "using transcripts in Read 2 for ${technology}"
            to=`echo $to | sed -e "s/_R3_/_I2_/g"`
        fi
        
        if [[ $verbose ]]; then
            echo "$to"
        fi
        crI2s+=($to)
        
        echo "handling $fq ..."
        if [[ ! -f $to ]] || [[ $convert  == true ]]; then
            if [[ $verbose ]]; then
                echo "cp -f $fq $to"
       	    fi
            cp -f $fq $to
        fi
    done
fi
##########



#####convert file format#####
echo "converting input files to confer cellranger format ..."
if [[ $convert == "false" ]]; then
    echo " input file format conversion skipped"
else
    echo " adjustment parameters:"
    echo "  barcodes: ${barcodeadjust}bps at its head"
    echo "  UMiIs: ${umiadjust}bps at its tail" 
    
    echo " making technology-specific modifications ..."
    #CEL-Seq2: swap barcode and UMI
    ## https://github.com/BUStools/bustools/issues/4
    if [[ "$technology" == "celseq2" ]]; then
        echo "  ...barcode and UMI swapped for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #swap UMI and barcode
            sed -E '2~2s/(.{6})(.{6})/\2\1/' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #inDrops: remove adapter (see links below for details)
    ## https://github.com/BUStools/bustools/issues/4
    ## https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html
    ## https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup
    #note that adapters do not have to be removed for the dual-indexed inDrops-v3
    if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
        echo "  ...remove adapter for ${technology}"
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
            echo"  ...barcode and UMI linker removed for ${technology}"
            sed -E '2~2s/^(.{8})(.{8}).{4}(.{6})/\1\2\3/g' > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #QuartzSeq: remove adapter
    if [[ "$technology" == "quartz-seq2-384" ]]; then
        for convFile in "${convFiles[@]}"; do
        echo "  ...remove adapter for ${technology}"
            #remove adapter if detected
            sed -E '
                /^TATAGAATTCGCGGCCGCTCGCGATAC(.{14})(.{8})/ {
                s/^TATAGAATTCGCGGCCGCTCGCGATAC(.{14})(.{8})/\1\2/g
                n
                n
                s/^.{27}(.{14})(.{8})/\1\2//g
                }'  $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    if [[ "$technology" == "quartz-seq2-1536" ]]; then
        echo "  ...remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapter if detected
            sed -E '
                /^TATAGAATTCGCGGCCGCTCGCGATAC(.{15})(.{8})/ {
                s/^TATAGAATTCGCGGCCGCTCGCGATAC(.{15})(.{8})/\1\2/g
                n
                n
                s/^.{27}(.{15})(.{8})/\1\2//g
                }'  $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #Sci-Seq: remove adapter and swap barcode and UMI
    if [[ "$technology" == "sciseq" ]]; then
        echo "  ...remove adapter for ${technology}"
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
            #swap barcode and UMI
            echo "  ...barcode and UMI swapped for ${technology}"
            sed -E '2~2s/(.{8})(.{10})/\2\1/' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile            
        done
    fi
    
    #SureCell: remove adapter and correct phase blocks
    ## https://github.com/Hoohm/dropSeqPipe/issues/42
    if [[ "$technology" == "surecell" ]]; then
        echo "  ...remove adapter and phase blocks for ${technology}"
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
    
    #SPLiT-Seq: correct phase blocks and swap barcode and UMI (if a whitelist and 18bp barcode can be supported)
    ## https://github.com/hms-dbmi/dropEst/issues/80
    ## https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup
    if [[ "$technology" == "splitseq" ]]; then
        echo "  ...remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers (swap barcode and UMI)
            sed -E '
                /^([ATCGA]{92})/ {
                s/^(.{10})(.{8}).{30}(.{8}).{30}(.{8})*/\2\3\4\1/g
                n
                n
                s/^(.{10})(.{8}).{30}(.{8}).{30}(.{8})*/\2\3\4\1/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        echo "  ...barcode and UMI swapped for ${technology}" #performed by \1 above
        done
    fi
    
    #SCRB-Seq: remove adapter
    ## https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html
    if [[ "$technology" == "scrbseq" ]]; then
        echo "  ...remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapters
            sed -E '
                /TCTTCCGATCT(.{6})(.{10})/ {
                s/TCTTCCGATCT(.{6})(.{10})/\1\2/g
                n
                n
                s/.{11}(.{6})(.{10})/\1\2/g
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
echo "cellranger count --id=$id\\
        --fastqs=$crIN\\
        --lanes=$LANE\\
        --r1-length=$totallength\\
        --chemistry=$chemistry\\
        --transcriptome=$reference\\
        --sample=$SAMPLE\\
        $d\\
        $n\\
        $j\\
        $l\\
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



#####process output#####
cloupefile=${SDIR}/${id}/outs/cloupe.cloupe
if [[ $technology != "10x" ]]; then
    #should not be necessary if setup is run correctly (will be omitted from Martian output)
    #this is kept to comply with the 10x End User License Agreement
    ###do not remove this code###
    if [[ -f $cloupefile ]]; then
        echo "Removing file ${cloupefile}"
        rm -rf $cloupefile
    fi
    echo "***Notice: Cloupe file cannot be computed for $technology"
    echo "           Cloupe files generated by this pipeline are corrupt"
    echo "           and cannot be read by the 10x Genomics Loupe Browser."
    echo "           We do not provide support for Cloupe files as this"
    echo "           requires software from 10x Genomics subject to their"
    echo "           End User License Agreement."
    echo "           Cloupe files are disabled in compliance with this."
fi
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
if [[ ${barcodefile} != "default:10x" ]]; then
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

