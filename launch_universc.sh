#!/bin/bash

#    UniverSC
#    Copyright (C) 2019  Tom Kelly; Kai Battenberg 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.


install=false

######UniverSC version#####
universcversion="1.2.7"
##########



#####locate cellranger and get cellranger version#####
cellrangerpath=$(which cellranger) #location of Cell Ranger
if [[ -z $cellrangerpath ]]; then
    echo "cellranger command is not found."
    exit 1
fi
## detects version by attempting to parse "cellranger --version" (older versions do not support this by print the help header if invalid args given)
cellrangerversion=$(cellranger --version 2>/dev/null | head -n 2 | tail -n 1 | rev | cut -f 1 -d" " | rev | cut -f 2 -d'(' | cut -f 1 -d')')
if [[ $verbose ]]; then
    echo $cellrangerversion
fi
##########


#####locate launch_universc.sh, barcode sources, and other tools######
echo "script running in $(realpath $0)..."
SOURCE=${BASH_SOURCE[0]}
# check if script is started via SLURM or bash (job ID defined in Slurm job)
if [[ -n $SLURM_JOB_ID ]];  then
    if [[ $(echo $SLURM_JOB_ID | wc -w) -ge 1 ]]; then
        # check the original location through scontrol and $SLURM_JOB_ID
        SOURCE=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
    fi
fi
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
echo "... script called from ${SDIR}"
RDIR="$( dirname "$SOURCE" )"

if [[ $RDIR != $SDIR ]]; then
    echo "'$RDIR' resolves to '$SDIR'"
fi
echo "Running launch_universc.sh in '$SDIR'"

TOOLS=${SDIR}/sub
ADDMOCKUMI=${TOOLS}/AddMockUMI.pl
BARCODERECOVER=${TOOLS}/RecoverBarcodes.pl
MAKEINDROPBARCODES=${TOOLS}/MakeIndropBarcodes.pl
FILTERSMARTSEQREADUMI=${TOOLS}/FilterSmartSeqReadUMI.pl
CONCATENATEBARCODES=${TOOLS}/ConcatenateDualIndexBarcodes.pl
PERCELLSTATS=${TOOLS}/ExtractBasicStats.pl
##########



#####define set options#####
# set defaults to 10x if missing (enables unit testing without warnings)
if [[ -z ${lastcall_b} ]]; then
    lastcall_b=16
fi
if [[ -z ${lastcall_u} ]]; then
    lastcall_u=10
fi
if [[ -z ${lastcall_p} ]]; then
    lastcall_p="default:10x"
fi
lockfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.lock #path for .lock file
lastcallfile=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.last_called #path for .last_called
lastcall=`[[ -e $lastcallfile ]] && cat $lastcallfile || echo ""`
lastcall_b=`echo ${lastcall} | cut -f 1 -d' '`
lastcall_u=`echo ${lastcall} | cut -f 2 -d' '`
lastcall_p=`echo ${lastcall} | cut -f 3 -d' '`
barcodedir=${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes #folder within cellranger with the whitelist barcodes
if [[ $(echo "${cellrangerversion} 4.0.0" | tr " " "\n" | sort -V | tail -n 1) == ${cellrangerversion} ]]; then
    barcodedir=$(dirname $cellrangerpath)/lib/python/cellranger/barcodes
    mkdir -p ${cellrangerpath}-cs/${cellrangerversion}
    ln -sf $(dirname $cellrangerpath)/mro/rna ${cellrangerpath}-cs/${cellrangerversion}/mro
    ln -sf $(dirname $cellrangerpath)/lib ${cellrangerpath}-cs/${cellrangerversion}/lib
fi
barcodefile=""
crIN=input4cellranger #name of the directory with all FASTQ and index files given to Cell Ranger
whitelistdir=${SDIR}/whitelists #path to whitelists
whitelistfile="outs/whitelist.txt" #name of the whitelist file added to the cellranger output
percellfile="outs/basic_stats.txt" #name of the file with the basic statistics of the run added to the cellranger output
##########



#####checking if Universc and Cell Ranger are writable#####
# user configuration
host=$(hostname) # returns host for local run and container ID for docker container
user=$(whoami 2> /dev/null || id -ur) # returns username if defined and user ID if not

#Cell Ranger
if ! [[ -w "$barcodedir" ]]; then
    echo "Error: Trying to run Cell Ranger installed at ${cellrangerpath}"
    echo "launch_universc.sh can only run with Cell Ranger installed locally"
    echo "Install Cell Ranger in a directory with write permissions such as /home/${user}/local and export to the PATH"
    echo "The following versions of Cell Ranger are found:"
    echo " `whereis cellranger`"
    exit 1
fi

#convert
if ! [[ -w "$SDIR" ]]; then
    echo "Error: Trying to run launch_universc.sh installed at $SDIR"
    echo "$SDIR must be writable to run launch_universc.sh"
    echo "Install launch_universc.sh in a directory with write permissions such as /home/${user}/local and export to the PATH"
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
    SHELL=$(readlink -f /proc/$$/exe | cut -d'/' -f 3)
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
if [[ $SHELL != "" ]]; then
    SHELL=" $SHELL"
fi
copyright="
    UniverSC Copyright (C) 2019 Tom Kelly; Kai Battenberg

"
notice="
    This program comes with ABSOLUTELY NO WARRANTY; for details type 'cat LICENSE'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type 'cat LICENSE' for details.

"
disclaimer="
    Cell Ranger is called as third-party dependency and is not maintained
    by this project. Please ensure you comply with the End User License
    Agreement for all software installed where applicable; for details type 'cat README.md'.

"
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
  -R1, --read1 FILE             Read 1 FASTQ file to pass to Cell Ranger (cell barcodes and umi)
  -R2, --read2 FILE             Read 2 FASTQ file to pass to Cell Ranger
  -I1, --index1 FILE            Index (I1) FASTQ file to pass to Cell Ranger (OPTIONAL)
  -I2, --index2 FILE            Index (I2) FASTQ file to pass to Cell Ranger (OPTIONAL and EXPERIMENTAL)
  -f,  --file NAME              Path and the name of FASTQ files to pass to Cell Ranger (prefix before R1 or R2)
                                  e.g. /path/to/files/Example_S1_L001
    
  -i,  --id ID                  A unique run id, used to name output folder
  -d,  --description TEXT       Sample description to embed in output files.
  -r,  --reference DIR          Path of directory containing 10x-compatible reference.
                                Available here the human genome and various model species:
                                https://genomec.gsc.riken.jp/gerg/UniverSC/Premade_references/

  -t,  --technology PLATFORM    Name of technology used to generate data.
                                Supported technologies:
                                  10x Genomics (version 2 or 3 automatically detected): 10x, chromium
                                  10x Genomics version 1 (14 bp barcode, 10 bp UMI): 10x-v1, chromium-v1
                                  10x Genomics version 2 (16 bp barcode, 10 bp UMI): 10x-v2, chromium-v2
                                  10x Genomics version 3 (16 bp barcode, 12 bp UMI): 10x-v3, chromium-v3
                                  Aligent Bravo B (16 bp barcode, No UMI): aligent, bravo
                                  BD Rhapsody v1 (27 bp barcode, 8 bp UMI): bd-rhapsody
                                  BD Rhapsody v2 enhanced beads (27 bp barcode, 8 bp UMI): bd-rhapsody-v2
                                  C1 Fluidigm (16 bp barcode, No UMI): c1, fluidgm-c1
                                  C1 CAGE (16 bp, No UMI): c1-cage
                                  C1 RamDA-Seq (16 bp, No UMI): c1-ramda-seq
                                  CEL-Seq (8 bp barcode, 4 bp UMI): celseq
                                  CEL-Seq2 (6 bp UMI, 6 bp barcode): celseq2
                                  Drop-Seq (12 bp barcode, 8 bp UMI): dropseq, nadia
                                  ICELL8 3\′ scRNA version 2 (11 bp barcode, No UMI): icell8-non-umi, icell8-v2
                                  ICELL8 3\′ scRNA version 3 (11 bp barcode, 14 bp UMI): icell8
                                  ICELL8 5\′ scRNA with TCR OR kit (10bp barcode, NO bp UMI): icell8-5-prime
                                  ICELL8 full-length scRNA with Smart-Seq (16 bp barcode, No UMI): icell8-full-length
                                  inDrops version 1 (19 bp barcode, 6 bp UMI): indrops-v1, 1cellbio-v1
                                  inDrops version 2 (19 bp barcode, 6 bp UMI): indrops-v2, 1cellbio-v2
                                  inDrops version 3 (16 bp barcode, 6 bp UMI): indrops-v3, 1cellbio-v3
                                  Nadia (12 bp barcode, 8 bp UMI): nadia, dropseq
                                  MARS-Seq (6 bp barcode, 10 bp UMI): marsseq, marsseq-v1
                                  MARS-Seq2 (7 bp barcode, 8 bp UMI): marsseq2, marsseq-v2
                                  Microwell-Seq (18 bp barcode, 6 bp UMI): microwell
                                  PIP-Seq version 0 (24 bp barcode, 8 bp UMI): pip-seq-v0
                                  PIP-Seq version 1 (16 bp barcode, 6 bp UMI): pip-seq-v1
                                  PIP-Seq version 2 (24 bp barcode, 12 bp UMI): pip-seq-v2
                                  PIP-Seq version 3 (28 bp barcode, 12 bp UMI): pip-seq-v3, fluent-bio-v3
                                  PIP-Seq version 4 (28 bp barcode, 12 bp UMI): pip-seq-v4, fluent-bio-v4
                                  QuartzSeq (6 bp index, no UMI): quartz-seq
                                  Quartz-Seq2 (14 bp barcode, 8 bp UMI): quartzseq2-384
                                  Quartz-Seq2 (15 bp barcode, 8 bp UMI): quartzseq2-1536
                                  RamDA-Seq (16 bp barcode, no UMI): ramda-seq
                                  SCI-Seq 2-level indexing (30 bp barcode, 8 bp UMI): sciseq2
                                  SCI-Seq 3-level indexing (40 bp barcode, 8 bp UMI): sciseq3
				  SCIFI-Seq (27 bp barcode, 8 bp UMI): scifiseq
                                  SCRB-Seq (6 bp barcode, 10 bp UMI): scrbseq, mcscrbseq
                                  SeqWell (12 bp barcode, 8 bp UMI): plexwell, seqwell, seqwells3
                                  Smart-seq (16 bp barcode, No UMI): smartseq
                                  Smart-seq2 (16 bp barcode, No UMI): smartseq2
                                  Smart-seq2-UMI, Smart-seq3 (16 bp barcode, 8 bp UMI): smartseq3
                                  SPLiT-Seq (8 bp UMI, 24 bp barcode): splitseq
                                  SPLiT-Seq v2.1 (10 bp UMI, 24 bp barcode): splitseq2
                                  STRT-Seq (6 bp barcode, no UMI): strt-seq
                                  STRT-Seq-C1 (8 bp barcode, 5 bp UMI): strt-seq-c1
                                  STRT-Seq-2i (13 bp barcode, 6 bp UMI): strt-seq-2i
                                  SureCell (18 bp barcode, 8 bp UMI): surecell, ddseq, biorad
                                  VASA-plate (6 bp UMI, 6 bp barcode): vasa-plate
                                  VASA-drop (6 bp UMI, 16 bp barcode): vasa-drop
                                Custom inputs are also supported by giving the name "custom" and length of barcode and UMI separated by "_"
                                  e.g. Custom (16 bp barcode, 10 bp UMI): custom_16_10
  
  -b,  --barcodefile FILE       Custom barcode list in plain text (with each line containing a barcode)
  
  -c,  --chemistry CHEM         Assay configuration, autodetection is not possible for converted files: 'SC3Pv2' (default), 'SC5P-PE', 'SC5P-R1', 'SC5P-R2', 'threeprime', or 'fiveprime'
                                    5\′ scRNA-Seq ('SC5P-PE') is available only for 10x Genomics, ICELL8, SmartSeq, and STRT-Seq technologies.
                                    Setting 'SC3Pv1' for 10x version 1 chemistry is recommended.
                                    All other technologies default to 3\′ scRNA-Seq parameters. Only 10x Genomics, ICELL8, and SmartSeq2 allow choosing which to use.
                                    For SmartSeq2 this parameter detemines using full-length sequences or 5\′ ends with internal reads removed.
  
  -n,  --force-cells NUM        Force pipeline to use this number of cells, bypassing the cell detection algorithm.
  -j,  --jobmode MODE           Job manager to use. Valid options: 'local' (default), 'sge', 'lsf', or a .template file
       --localcores NUM         Set max cores the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --localmem NUM           Set max GB the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --mempercore NUM         Set max GB each job may use at one time.
                                    Only applies in cluster jobmodes.
  
  -p,  --per-cell-data          Generates a file with basic run statistics along with per-cell data
  
       --non-umi or --read-only Force counting reads by adding or replacing UMI with a mock sequence.
                                Default for: Quartz-Seq, RamDA-Seq, Smart-Seq, Smart-Seq2, STRT-Seq, ICELL8 3-prime version2
  
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

echo $copyright

echo $notice

echo $disclaimer

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
bam=""
secondary=""
library=""
introns=""
noexit=""
nopreflight=""
nonUMI=false
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
        --non-umi|--read-only)
            nonUMI=true
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
        --include-introns)
            introns="include-introns"
            next=false
            shift
            ;;
        --no-bam)
            bam="--no-bam"
            next=false
            shift
            ;;
        --nosecondary)
            secondary="--nosecondary"
            next=false
            shift
            ;;
        --nolibraries)
            library="--nolibraries"
            next=false
            shift
            ;;
        --noexit)
            noexit="--noexit"
            next=false
            shift
            ;;
        --nopreflight)
            nopreflight="--nopreflight"
            next=false
            shift
            ;;
        -h|--help)
            echo "$help"
            exit 0
            ;;
        -v|--version)
            echo "launch_universc.sh version ${universcversion}"
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
        ncells=`wc -l $barcodefile | cut -f 1 -d' '`
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
elif [[ "$technology" == "10x-v1" ]] || [[ "$technology" == "chromium-v1" ]]; then
     technology="10x-v1"
elif [[ "$technology" == "10x-v2" ]] || [[ "$technology" == "chromium-v2" ]]; then
    technology="10x-v2"
elif [[ "$technology" == "10x-v3" ]] || [[ "$technology" == "chromium-v3" ]]; then
    technology="10x-v3"
elif [[ "$technology" == "bd-rhapsody" ]] || [[ "$technology" == "bd" ]] || [[ "$technology" == "rhapsody" ]] || [[ "$technology" == "bdrhapsody" ]]; then
    technology="bd-rhapsody"
elif [[ "$technology" == "bd-rhapsody-v2" ]] || [[ "$technology" == "bd-v2" ]] || [[ "$technology" == "rhapsody-v2" ]] || [[ "$technology" == "bdrhapsody-v2" ]] || [[ "$technology" == "bd-rhapsodyv2" ]] || [[ "$technology" == "bdv2" ]] || [[ "$technology" == "rhapsodyv2" ]] || [[ "$technology" == "bdrhapsodyv2" ]] || [[ "$technology" == "bd-rhapsody-2022" ]] || [[ "$technology" == "bd-2022" ]] || [[ "$technology" == "rhapsody-2022" ]] || [[ "$technology" == "bdrhapsody-2022" ]] || [[ "$technology" == "bd-rhapsody2022" ]] || [[ "$technology" == "bd2022" ]] || [[ "$technology" == "rhapsody2022" ]] || [[ "$technology" == "bdrhapsody2022" ]]; then
    technology="bd-rhapsody-v2"
elif [[ "$technology" == "aligent" ]] || [[ "$technology" == "bravo" ]] || [[ "$technology" == "aligent-bravo" ]] || [[ "$technology" == "bravo-b" ]] || [[ "$technology" == "aligent-bravo-b" ]]; then
    techology="bravo"
    nonUMI=false
elif [[ "$technology" == "c1" ]] || [[ "$technology" == "c1-fluidigm" ]] || [[ "$technology" == "fluidigm" ]] || [[ "$technology" == "fluidigm-c1" ]]|| [[ "$technology" == "fluidigmc1" ]] ||  [[ "$technology" == "c1-rna-seq" ]]|| [[ "$technology" == "c1-mrna-seq" ]] || [[ "$technology" == "c1-rnaseq" ]]|| [[ "$technology" == "c1-scrna" ]]; then
    technology="fluidigm-c1"
    nonUMI=true
elif [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "c1cage" ]] || [[ "$technology" == "cage-c1" ]] || [[ "$technology" == "cagec1" ]]; then
    technology="c1-cage"
    nonUMI=true
elif [[ "$technology" == "celseq" ]] || [[ "$technology" == "cel-seq" ]]; then
    technology="celseq"
elif [[ "$technology" == "celseq2" ]] || [[ "$technology" == "cel-seq2" ]]; then
    technology="celseq2"
elif [[ "$technology" == "nadia" ]] || [[ "$technology" == "dropseq" ]] || [[ "$technology" == "drop-seq" ]]; then
    technology="nadia"
elif [[ "$technology" == "icell8" ]] || [[ "$technology" == "icell-8" ]] ||  [[ "$technology" == "icell8-v3" ]] || [[ "$technology" == "icell8v3" ]]; then
    technology="icell8"
    #set as "icell8-5-prime" if called by chemistry
    if [[ "$chemistry" == "SC5P-"* ]] || [[ "$chemistry" == "fiveprime" ]]; then
        technology="icell8-icell8-5-prime"
        nonUMI=true
    fi
elif [[ "$technology" == "icell8-non-umi" ]] || [[ "$technology" == "icell8-nonumi" ]] || [[ "$technology" == "icell8-v2" ]] || [[ "$technology" == "icell8v2" ]]; then
    technology="icell8"
    nonUMI=true
elif [[ "$technology" == "icell8-five-prime" ]] || [[ "$technology" == "icell8fiveprime" ]] || [[ "$technology" == "icell8-5p" ]] || [[ "$technology" == "icell85p" ]]; then
    technology="icell8-5-prime"
    nonUMI=true
    if [[ -z ${chemistry} ]] || [[ ${chemistry} == "SC3Pv"* ]]; then
        if [[ $verbose ]]; then
            echo "setting chemistry for 5' scRNA"
        fi
        chemistry="SC5P-PE"
    fi
elif [[ "$technology" == "icell8-full" ]] || [[ "$technology" == "icell8-fl" ]] || [[ "$technology" == "icell8fl" ]] || [[ "$technology" == "smartseq-icell8" ]] || [[ "$technology" == "smart-seq-icell8" ]] || [[ "$technology" == "icell8-smart-seq" ]] || [[ "$technology" == "icell8smartseq" ]]; then
    technology="icell8-full-length"
    nonUMI=true
    if [[ -z ${chemistry} ]] || [[ ${chemistry} == "SC3Pv"* ]]; then
        if [[ $verbose ]]; then
            echo "setting chemistry for 5' scRNA"
        fi
        chemistry="SC5P-PE"
    fi
elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrops-v1" ]] || [[ "$technology" == "indropv1" ]] || [[ "$technology" == "indropsv1" ]] || [[ "$technology" == "indrop1" ]] || [[ "$technology" == "indrops1" ]] || [[ "$technology" == "1cellbio-v1" ]] || [[ "$technology" == "1cellbiov1" ]]; then
    technology="indrop-v1"
elif [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrops-v2" ]] || [[ "$technology" == "indropv2" ]] || [[ "$technology" == "indropsv2" ]] || [[ "$technology" == "indrop2" ]] || [[ "$technology" == "indrops2" ]] ||[[ "$technology" == "1cellbio-v2" ]] || [[ "$technology" == "1cellbiov2" ]]; then
    technology="indrop-v2"
elif [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "indrops-v3" ]] || [[ "$technology" == "indropv3" ]] || [[ "$technology" == "indropsv3" ]] || [[ "$technology" == "indrop3" ]] || [[ "$technology" == "indrops3" ]] || [[ "$technology" == "1cellbio-v3" ]] || [[ "$technology" == "1cellbiov3" ]]; then
    technology="indrop-v3"
elif [[ "$technology" == "marsseq" ]] || [[ "$technology" == "mars-seq" ]] || [[ "$technology" == "marsseq-v1" ]] || [[ "$technology" == "mars-seq-v1" ]] || [[ "$technology" == "marsseqv1" ]] || [[ "$technology" == "mars-seqv1" ]]; then
    technology="marsseq-v1"
elif [[ "$technology" == "marsseq2" ]] || [[ "$technology" == "mars-seq2" ]] || [[ "$technology" == "marsseq-v2" ]] || [[ "$technology" == "mars-seq-v2" ]] || [[ "$technology" == "marsseqv2" ]] || [[ "$technology" == "mars-seqv2" ]]; then
    technology="marsseq-v2"
elif [[ "$technology" == "microwell-seq" ]] || [[ "$technology" == "micro-well" ]] || [[ "$technology" == "microwell" ]] || [[ "$technology" == "microwellseq" ]]; then
     technology="microwellseq"
elif [[ "$technology" == "pip-seq-v0" ]] || [[ "$technology" == "pip-seq-V0" ]] || [[ "$technology" == "pipseq-v0" ]] || [[ "$technology" == "pip-seqv0" ]] || [[ "$technology" == "pipseqv0" ]] || [[ "$technology" == "delley" ]] || [[ "$technology" == "delley2021" ]] || [[ "$technology" == "delley-2021" ]]; then
     technology="pip-seq-v0"
elif [[ "$technology" == "pip-seq" ]] || [[ "$technology" == "pip-seq-v1" ]] || [[ "$technology" == "pip-seq" ]] || [[ "$technology" == "pip-seq-v1" ]] || [[ "$technology" == "pipseq" ]] || [[ "$technology" == "pipseq-v1" ]] || [[ "$technology" == "pip-seqv1" ]] || [[ "$technology" == "pipseqv1" ]]; then
     technology="pip-seq-v1"
elif [[ "$technology" == "pip-seq-v2" ]] || [[ "$technology" == "pip-seq-V2" ]] || [[ "$technology" == "pipseq-v2" ]] || [[ "$technology" == "pip-seqv2" ]] || [[ "$technology" == "pipseqv2" ]]; then
     technology="pip-seq-v2"
elif [[ "$technology" == "pip-seq-v3" ]] || [[ "$technology" == "pip-seq-V3" ]] || [[ "$technology" == "pipseq-v3" ]] || [[ "$technology" == "pip-seqv3" ]] || [[ "$technology" == "pipseqv3" ]] || [[ "$technology" == "fluent-bio-v3" ]] || [[ "$technology" == "fluentbio-v3" ]] || [[ "$technology" == "fluent-biov3" ]] || [[ "$technology" == "fluentbiov3" ]]; then
     technology="pip-seq-v3"
elif [[ "$technology" == "pip-seq-v4" ]] || [[ "$technology" == "pip-seq-V4" ]] || [[ "$technology" == "pipseq-v4" ]] || [[ "$technology" == "pip-seqv4" ]] || [[ "$technology" == "pipseqv4" ]] || [[ "$technology" == "fluent-bio-v4" ]] || [[ "$technology" == "fluentbio-v4" ]] || [[ "$technology" == "fluent-biov4" ]] || [[ "$technology" == "fluentbiov4" ]]; then
     technology="pip-seq-v4"
elif [[ "$technology" == "quartz-seq" ]] || [[ "$technology" == "quartzseq" ]] || [[ "$technology" == "quartz-seq1" ]]; then
     technology="quartz-seq"
     nonUMI=true
elif [[ "$technology" == "quartz-seq2-384" ]] || [[ "$technology" == "quartzseq2-384" ]] || [[ "$technology" == "quartz-seq2-v3.1" ]] || [[ "$technology" == "quartzseq2-v3.1" ]] || [[ "$technology" == "quartzseq2v3.1" ]]; then
    technology="quartz-seq2-384"
elif [[ "$technology" == "quartz-seq2-1536" ]] || [[ "$technology" == "quartzseq2-1536" ]] || [[ "$technology" == "quartz-seq2-v3.2" ]] || [[ "$technology" == "quartzseq2-v3.2" ]] || [[ "$technology" == "quartzseq2v3.2" ]]; then
    technology="quartz-seq2-1536"
elif [[ "$technology" == "rambda-seq" ]] || [[ "$technology" == "lamda-seq" ]] || [[ "$technology" == "lambda-seq" ]] || [[ "$technology" == "ramdaseq" ]] || [[ "$technology" == "ram-da-seq" ]] || [[ "$technology" == "ramda-seq" ]]; then
     technology="ramda-seq"
     nonUMI=true
elif [[ "$technology" == "rambda-seq-c1" ]] || [[ "$technology" == "lamda-seq-c1" ]] || [[ "$technology" == "lambda-seq-c1" ]] || [[ "$technology" == "ramdaseqc1" ]] || [[ "$technology" == "ram-da-seq-c1" ]] || [[ "$technology" == "ramda-seq-c1" ]] || [[ "$technology" == "ramda-seqc1" ]] || [[ "$technology" == "c1-rambda-seq" ]] || [[ "$technology" == "c1-lamda-seq" ]] || [[ "$technology" == "c1-lambda-seq" ]] || [[ "$technology" == "c1-ramdaseq" ]] || [[ "$technology" == "c1-ram-da-seq" ]] || [[ "$technology" == "c1ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]]; then
     technology="c1-ramda-seq"
     nonUMI=true
elif [[ "$technology" == "sciseq" ]] || [[ "$technology" == "sci-seq" ]] || [[ "$technology" == "sci-rna-seq" ]]; then
    technology="sciseq3"
elif [[ "$technology" == "sciseq2" ]] || [[ "$technology" == "sci-seq2" ]]; then
     technology="sciseq2"
elif [[ "$technology" == "sciseq3" ]] || [[ "$technology" == "sci-seq3" ]] || [[ "$technology" == "sci-rna-seq3" ]] || [[ "$technology" == "sci-rna-seq-3" ]] ; then
     technology="sciseq3"
elif [[ "$technology" == "scifiseq" ]] || [[ "$technology" == "scifi-seq" ]] || [[ "$technology" == "sci-fi-seq" ]] || [[ "$technology" == "scifi-rna-seq" ]]; then
      technology="scifiseq"
elif [[ "$technology" == "scrbseq" ]] || [[ "$technology" == "scrb-seq" ]] || [[ "$technology" == "mcscrbseq" ]] || [[ "$technology" == "mcscrb-seq" ]]; then
    technology="scrbseq"
elif [[ "$technology" == "plexwell" ]] || [[ "$technology" == "plex-well" ]] ||  [[ "$technology" == "seqwell" ]] || [[ "$technology" == "seq-well" ]] || [[ "$technology" == "seqwells3" ]] || [[ "$technology" == "seq-well-s3" ]]; then
    technology="seqwell"
elif [[ "$technology" == "smartseq" ]] || [[ "$technology" == "smart-seq" ]] || [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smart-seq2" ]]; then
    technology="smartseq2"
    nonUMI=true
elif [[ "$technology" == "smartseq2-umi" ]] || [[ "$technology" == "smart-seq2-umi" ]]; then
    technology="smartseq2-umi"
    nonUMI=false
elif [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "smart-seq3" ]]; then
    technology="smartseq3"
    nonUMI=false
elif [[ "$technology" == "splitseq" ]] || [[ "$technology" == "split-seq" ]]; then
    technology="splitseq"
elif [[ "$technology" == "splitseq2" ]] || [[ "$technology" == "split-seq2" ]] || [[ "$technology" == "splitseq-v2" ]] || [[ "$technology" == "split-seq-v2" ]] || [[ "$technology" == "parse" ]] || [[ "$technology" == "parsebio" ]] || [[ "$technology" == "parse-bio" ]] || [[ "$technology" == "evercode" ]]; then
    technology="splitseq2"
elif [[ "$technology" == "strt-seq" ]] || [[ "$technology" == "strt" ]] || [[ "$technology" == "strtseq" ]]; then
    technology="strt-seq"
    nonUMI=true
    if [[ -z ${chemistry} ]] || [[ ${chemistry} == "SC3Pv"* ]]; then
        if [[ $verbose ]]; then
            echo "setting chemistry for 5\' scRNA"
        fi
        chemistry="SC5P-R1"
    fi
elif [[ "$technology" == "strt-seq-c1" ]] || [[ "$technology" == "strt-seqc1" ]] || [[ "$technology" == "strtseqc1" ]] || [[ "$technology" == "strtseq-c1" ]]; then
     technology="strt-seq-c1"
elif [[ "$technology" == "strt-seq-2i" ]] || [[ "$technology" == "strt-seq2i" ]] || [[ "$technology" == "strtseq2i" ]] || [[ "$technology" == "strtseq-2i" ]]; then
     technology="strt-seq-2i"
elif [[ "$technology" == "strt-seq-2018" ]] || [[ "$technology" == "strt-seqc2018" ]] || [[ "$technology" == "strtseq2018" ]] || [[ "$technology" == "strtseq-2018" ]] || [[ "$technology" == "strt-seq-v3" ]]; then
    technology="strt-seq-2018"
    if [[ -z ${chemistry} ]] || [[ ${chemistry} == "SC5P"* ]]; then
        if [[ $verbose ]]; then
            echo "setting chemistry for 3\' scRNA PE"
        fi
        chemistry="SC3Pv2"
     fi
elif [[ "$technology" == "surecell" ]] || [[ "$technology" == "surecellseq" ]] || [[ "$technology" == "surecell-seq" ]] || [[ "$technology" == "ddseq" ]] || [[ "$technology" == "dd-seq" ]] || [[ "$technology" == "bioraad" ]]; then
    technology="surecell"
elif [[ "$technology" == "vasa-seq-plate" ]] || [[ "$technology" == "vasa-seq" ]] || [[ "$technology" == "vasa-plate" ]] || [[ "$technology" == "vasaplate" ]] || [[ "$technology" == "vasa" ]]; then
    technology="vasa-plate"
elif [[ "$technology" == "vasa-seq-drop" ]] || [[ "$technology" == "vasa-drop-seq" ]] || [[ "$technology" == "vasa-drop" ]] || [[ "$technology" == "vasadrop" ]] || [[ "$technology" == "vasa-droplet" ]]; then
    technology="vasa-drop"
elif [[ "$technology" == "custom"* ]]; then
    fields=$((`echo $technology | grep -o "_" | wc -l` + 1))
    if [[ $fields -ne 3 ]]; then
        echo "Error: custom input must have exactly 3 fields separated by '_', e.g. custom_10_16"
    else
        b=`echo $technology | cut -f 2 -d'_'`
        u=`echo $technology | cut -f 3 -d'_'`
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
    echo "***WARNING: ${technology}-v3 should only be used for kits that have valid UMIs***"
fi
if [[ "$technology" == "smartseq2-umi" ]] || [[ "$technology" == "smartseq3" ]]; then
    echo "***WARNING: ${technology} should only be used for kits that have UMIs***"
    echo "... UMI reads will be filtered using a tag sequence which will be subsequently removed"
    echo "... barcodes will be derived from dual indexes"
fi
if [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
    echo "***Note: Make sure that samples are demultiplexed prior to running launch_universc.sh***"
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
elif [[ "$technology" == "10x-v1" ]]; then
    barcodelength=14
    barcode_default=14
    umilength=10
    minlength=14
    chemistry="SC3Pv1"
    technology="10x"
elif [[ "$technology" == "10x-v2" ]]; then
    barcodelength=16
    umilength=10
    minlength=16
    technology="10x"
elif [[ "$technology" == "10x-v3" ]]; then
    barcodelength=16
    umilength=12
    minlength=16
    chemistry="SC3Pv3"
    technology="10x"
elif [[ "$technology" == "bd-rhapsody" ]] || [[ "$technology" == "bd-rhapsody-v2" ]]; then
    barcodelength=27
    umilength=8
    minlength=27
elif [[ "$technology" == "fluidigm-c1" ]]; then
     barcodelength=16
     umilength=0
     minlength=16
elif [[ "$technology" == "c1-cage" ]]; then
     barcodelength=16
     umilength=0
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
    if [[ $nonUMI == "true" ]]; then
       umilength=0
    else
       umilength=14
    fi
    minlength=11
elif [[ "$technology" == "icell8-5-prime" ]]; then
    barcodelength=10
    if [[ $nonUMI == "true" ]]; then
        umilength=0
    else
        umilength=8
    fi
    minlength=10
elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
    barcodelength=19
    umilength=6
    minlength=16
elif [[ "$technology" == "indrop-v3" ]]; then
    barcodelength=16
    umilength=6
    minlength=16
elif [[ "$technology" == "marsseq-v1" ]]; then
    barcodelength=6
    umilength=10
    minlength=6
elif [[ "$technology" == "microwellseq" ]]; then
    barcodelength=18
    umilength=6
    minlength=18
elif [[ "$technology" == "marsseq-v2" ]]; then
    barcodelength=7
    umilength=8
    minlength=7
elif [[ "$technology" == "pip-seq-v0" ]]; then
    barcodelength=27
    umilength=8
    minlength=24
elif [[ "$technology" == "pip-seq-v1" ]]; then
    barcodelength=19
    umilength=6
    minlength=16
elif [[ "$technology" == "pip-seq-v2" ]]; then
    barcodelength=24
    umilength=12
    minlength=24
elif [[ "$technology" == "pip-seq-v3" ]]; then
    barcodelength=28
    umilength=12
    minlength=28
elif [[ "$technology" == "pip-seq-v4" ]]; then
    barcodelength=31
    umilength=12
    minlength=28
elif [[ "$technology" == "quartz-seq2-384" ]]; then
    barcodelength=14
    umilength=8
    minlength=14
elif [[ "$technology" == "quartz-seq2-1536" ]]; then
    barcodelength=15
    umilength=8
    minlength=15
elif [[ "$technology" == "sciseq2" ]]; then
    barcodelength=30
    umilength=8
    minlength=30
elif [[ "$technology" == "sciseq3" ]]; then
     barcodelength=40
     umilength=8
     minlength=40
elif [[ "$technology" == "scifiseq" ]]; then
     barcodelength=27
     umilength=8
     minlength=27
elif [[ "$technology" == "scrbseq" ]]; then
    barcodelength=6 
    umilength=10
    minlength=6
elif [[ "$technology" == "seqwell" ]]; then
    barcodelength=8
    umilength=12
    minlength=8
elif [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "icell8-full-length" ]] || [[ "$technology" == "fluidigm-c1" ]] || [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "quartz-seq" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]] || [[ "$technology" == "bravo" ]]; then
    barcodelength=16
    if [[ $nonUMI == "true" ]]; then
       umilength=0
    else
       umilength=8
    fi
    minlength=16
elif [[ "$technology" == "splitseq" ]]; then
     barcodelength=24
     umilength=8
     minlength=24
elif [[ "$technology" == "splitseq2" ]]; then
     barcodelength=24
     umilength=10
     minlength=24
elif [[ "$technology" == "strt-seq" ]]; then
    barcodelength=6
    umilength=0
    minlength=6
elif [[ "$technology" == "strt-seq-c1" ]]; then
    barcodelength=8
    umilength=5
    minlength=8
elif [[ "$technology" == "strt-seq-2i" ]]; then
    barcodelength=13
    umilength=6
    minlength=13
elif [[ "$technology" == "strt-seq-2018" ]]; then
    barcodelength=8
    umilength=8
    minlength=8
elif [[ "$technology" == "surecell" ]]; then
    barcodelength=18
    umilength=8
    minlength=18
elif [[ "$technology" == "vasa-plate" ]]; then
    barcodelength=6
    umilength=6
    minlength=6
elif [[ "$technology" == "vasa-drop" ]]; then
    barcodelength=16
    umilength=6
    minlength=16
elif [[ "$technology" == "custom"* ]]; then
    barcodelength=`echo $technology | cut -f 2 -d'_'`
    umilength=`echo $technology | cut -f 3 -d'_'`
    minlength=${barcodelength}
fi
if [[ $verbose ]]; then
    echo " barcode and UMI lengths set as ${barcodelength} and ${umilength} respectively"
fi
##########



#####Setting chemistry#####
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
    [[ "$chemistry" != "SC5P-R1" ]] && \
    [[ "$chemistry" != "SC5P-R2" ]] && \
    [[ "$chemistry" != "SC-FB" ]]; then
        echo "Error: option -c can be auto, threeprime, fiveprime, SC3Pv1, SC3Pv2, SC3Pv3, SC5P-PE, SC5P-R1, SC5P-R2, or SC-FB"
        exit 1
    fi
fi

#determine what chemistry is recommended
temp_chemistry="SC3Pv2"
if [[ $umilength -gt 10 ]]; then
    temp_chemistry="SC3Pv3"
fi
if [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "strt-seq"* ]]; then
    temp_chemistry="SC5P-R1"
fi
if [[ "$technology" == "10x" ]]; then
    temp_chemistry="auto"
fi

#set chemistry
if [[ -z ${chemistry} ]]; then
    chemistry=${temp_chemistry}
elif [[ "$chemistry" != "$temp_chemistry" ]]; then
    echo "***WARNING: chemistry is set to ${chemistry} where ${temp_chemistry} would have been chosen automatically. proceed with caution.***"
fi

#set default barcode and umi lengths
if [[ $minlength -gt 16 ]]; then
    barcode_default=$minlength
elif [[ "$chemistry" == "SC3Pv1" && ${technology} == "10x" ]]; then
    barcode_default=14
else
    barcode_default=16
fi
if [[ $umilength -gt 12 ]]; then
    umi_default=$umilength
else
    if [[ "$chemistry" == "SC3Pv3" ]]; then
        umi_default=12
    else
        umi_default=10
    fi
fi

totallength=`echo $((${barcode_default} + ${umi_default}))`

#adjustment lengths
if [[ $verbose ]]; then
   echo "barcode length: $barcodelength"
   echo "barcode default: $barcode_default"
fi
barcodeadjust=`echo $((${barcodelength} - ${barcode_default}))`
umiadjust=`echo $((${umilength} - ${umi_default}))`
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
    elif [[ ${#read2[@]} -eq 0 ]] && [[ "$chemistry" != "SC5P-R1" ]]; then
        echo "Error: option --read2 or --file is required"
        exit 1
    fi
fi

#check for presence of index files
if [[ $verbose ]]; then
    echo " checking option: --index1 and --index2"
fi


#index 1
r4_present="false"
i1_present="false"
if [[ $setup == "false" ]]; then
    if [[ ${#index1[@]} -ne ${#read1[@]} ]]; then
        if [[ ${#index1[@]} -gt 0 ]]; then
            echo " Error: number of index1 files is not matching the number of R1 and R2 files"
            echo " either give no index1 file, or give index1 file for each and every read1 file"
            exit 1
        else
            if [[ $verbose ]]; then
                echo " No index files given. Automatically detecting from R1 and R2 file names ..."
            fi
            r1_list=("${read1[@]}")
            r2_list=("${read2[@]}")
            i1_list=()
            r4_list=()
            for j in ${!r1_list[@]}; do
                read=${r1_list[$j]}
                R1_file=$read
                if [[ $read == *R1* ]]; then
                    R2_file=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
                    R4_file=$(echo $read | perl -pne 's/(.*)_R1/$1_R4/' )
                    I1_file=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
                    if [[ -f $R4_file ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R4_file})'*.gz' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R4_file})'*.fastq' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R4_file})'*.fq' | head -n 1) ]]; then
                        r4_present="true"
                        if [[ $verbose ]]; then
                            echo "  file $R4_file found, replacing $R1_file ..." 
                        fi
                        r1_read=$R4_file
                        r1_list[$j]=$r1_read
                        if [[ $verbose ]]; then
                            echo "  file $R1_file found, replacing $R2_file ..."
                        fi
                        r2_read=$R1_file
                        r2_list[$j]=$r2_read
                        if [[ $verbose ]]; then
                            echo "  file $R2_file found, replacing $I1_file ..."
                        fi
                        i1_read=$R2_file
                        i1_list[$j]=$i1_read
                    fi
                    if [[ -f $I1_file ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I1_file})'*.gz' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I1_file})'*.fastq' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I1_file})'*.fq' | head -n 1) ]]; then
                        i1_present="true"
                        if [[ $verbose ]]; then
                            echo "  file $I1_file found ..."
                        fi
                        i1_read=$I1_file
                        i1_list[$j]=$i1_read
                    fi
                fi
            done
            read1=("${r1_list[@]}")
            read2=("${r2_list[@]}")
            if [[ ${#i1_list[@]} -eq 0 ]]; then
                if [[ $verbose ]]; then
                    echo "No index files found"
                fi
            else
                index1=("${i1_list[@]}")
                if [[ $verbose ]]; then
                    echo "  index files found ${#index1[@]} I1s: ${index1[@]}"
                fi
            fi
        fi
    fi
fi

#index 2
r3_present="false"
i2_present="false"
if [[ $setup == "false" ]]; then
    #only check I2 for dual-indexed techniques
    if [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "sciseq2" ]] || [[ "$technology" == "sciseq3" ]] || [[ "$technology" == "scifiseq" ]] || [[ "$technology" == "smartseq"* ]] || [[  "$chemistry" == "SC3Pv1" && "$technology" == "10x" ]]; then
        if [[ ${#index2[@]} -ne ${#read1[@]} ]]; then
            if [[ ${#index2[@]} -gt 0 ]]; then
               echo " Error: number of index1 files is not matching the number of index2 files"
               echo " for $technology, either give no index files or give index1 and index2 for each and every read1 and read2 file"
               exit 1
            else
                if [[ $verbose ]]; then
                    echo " No index files given. Automatically detecting from R1 and R2 file names ..."
                fi
                r1_list=("${read1[@]}")
                r2_list=("${read2[@]}")
                r3_list=()
                i2_list=()
                for j in ${!r1_list[@]}; do
                    read=${r1_list[$j]}
                    R1_file=$read
                    if [[ $read == *R1* ]]; then
                        R2_file=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
                        R3_file=$(echo $read | perl -pne 's/(.*)_R1/$1_R3/' )
                        I2_file=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
                        if [[ -f $R3_file ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R3_file})'*.gz' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R3_file})'*.fastq' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${R3_file})'*.fq' | head -n 1) ]]; then
                            r3_present="true"
                            if [[ $verbose ]]; then
                                echo "  file $R3_file found, replacing $I2_file ..." 
                            fi
                            if [[ $technology == "10x" ]]; then
                                r3_read=$R3_file
                                r3_list[$j]=$r3_read
                            else
                                i2_read=$R3_file
                                i2_list[$j]=$r3_read
                            fi
                        fi
                        if [[ -f $I2_file ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I2_file})'*.gz' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I2_file})'*.fastq' | head -n 1) ]] || [[ -f $(find $(dirname ${read}) -name $(basename ${I2_file})'*.fq' | head -n 1) ]]; then
                            i2_present="true"
                            if [[ $verbose ]]; then
                                echo "  file $I2_file found ..."
                            fi
                            i2_read=$I2_file
                            i2_list[$j]=$i2_read
                        fi
                    fi
                done
                if [[ ${#i2_list[@]} -eq 0 ]]; then
                    if [[ $verbose ]]; then
                        echo "Dual index files not found"
                    fi
                else
                    index2=("${i2_list[@]}")
                    if [[ $verbose ]]; then
                        echo "  index files found ${#index2[@]} I2s: ${index2[@]}"
                    fi
                fi
                if [[ ${#r3_list[@]} -gt 0 ]]; then
                    read3=("${r3_list[@]}")
                    if [[ $verbose ]]; then
                        echo "  index files found ${#read3[@]} R3s: ${read3[@]}"
                    fi
                fi
            fi
        fi
    elif [[ ${#index2[@]} -gt 0 ]]; then
        echo " Error: $technology does not support dual index"
        echo " re-run without selecting any files for --index2"
        exit 1
    fi
fi

if [[ $verbose ]]; then
    echo " Input files post-curation 1"
    echo "  ${#read1[@]}files - read1s: ${read1[@]}"
    echo "  ${#read2[@]}files - read2s: ${read2[@]}"
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#read3[@]}files - R3s: ${read3[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#read4[@]}files - R4s: ${read4[@]}"
    fi
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]}files - I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]}files - I2s: ${index2[@]}"
    fi
    echo "  number of these files are as expected"
fi



keys=("R1" "R2")
if [[ $r3_present == "true" ]]; then
    keys=(${keys[@]} "R3")
    if [[ $technology == "10x" ]]; then
        read3=("${read3[@]}")
    else
        read3=("${index2[@]}")
    fi
fi
if [[ $r4_present == "true" ]]; then
    keys=(${keys[@]} "R4")
    read4=("${index1[@]}")
fi
if [[ $i1_present == "true" ]]; then
    keys=(${keys[@]} "I1")
fi
if [[ $i2_present == "true" ]]; then
    keys=(${keys[@]} "I2")
fi
##########



#####Input file curation 2: Check R1, R2, I1, and I2 files for their extensions#####
##allows incomplete file names and processing compressed files
if [[ $verbose ]]; then
    echo "keys: ${keys[@]}"
fi
for key in ${keys[@]}; do
    if [[ $verbose ]]; then
        echo "key: $key"
    fi
    readkey=$key
    list=()
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    elif [[ $readkey == "R3" ]]; then
        list=("${read3[@]}")
    elif [[ $readkey == "R4" ]]; then
        list=("${read4[@]}")
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
                echo "  found compressed file ..."
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
        elif [[ -f ${read}.gz ]] || [[ -h ${read}.gz ]]; then
            read=${read}.gz
            gunzip -f -k ${read}
            read=`echo $read | sed -e "s/\.gz//g"`
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
    elif [[ $readkey == "R3" ]]; then
        read3=("${list[@]}")
    elif [[ $readkey == "R4" ]]; then
        read4=("${list[@]}")
    elif [[ $readkey == "I1" ]]; then
        index1=("${list[@]}")
     elif [[ $readkey == "I2" ]]; then
        index2=("${list[@]}")
    fi
done

if [[ $verbose ]]; then
    echo " Input files post-curation 2"
    echo "  ${#read1[@]}files - read1s: ${read1[@]}"
    echo "  ${#read2[@]}files - read2s: ${read2[@]}"
    if [[ ${#read3[@]} -gt 0 ]]; then
        echo "  ${#read3[@]}files - R3s: ${read3[@]}"
    fi
    if [[ ${#read4[@]} -gt 0 ]]; then
        echo "  ${#read4[@]}files -  R4s: ${read4[@]}"
    fi
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]}files - I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]}files -  I2s: ${index2[@]}"
    fi
    echo "  files exist, with extentions compatible with launch_universc.sh"
fi
##########



#####Input file curation 3: renaming read1, read2, index1, and index2 file name if not compatible with the launch_universc.sh#####
for key in ${keys[@]}; do
    if [[ $verbose ]]; then
        echo "key: $key"
    fi
    readkey=$key
    list=""
    if [[ $readkey == "R1" ]]; then
        list=("${read1[@]}")
    elif [[ $readkey == "R2" ]]; then
        list=("${read2[@]}")
    elif [[ $readkey == "R3" ]]; then
        list=("${read3[@]}")
    elif [[ $readkey == "R4" ]]; then
        list=("${read4[@]}")
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
            fullpath=`readlink -f $read`
            if [[ $verbose ]]; then
                echo "***Warning: file $read not in current directory. Path to the file captured instead.***"
                echo " (file) $read"
                echo " (path) $fullpath"
            fi
            read=${fullpath}
        fi
        if [[ $read == *$readkey* ]]; then
            if [[ $verbose ]]; then
                echo " (file) $read contains $readkey"
            fi
        else
            if [[ $verbose ]]; then
                echo " (file) $read does not contain $readkey"
            fi
            readsuffix="${readkey: -1}"
            rename -f "s/_${readsuffix}\./_R${readsuffix}_001./" $read
            read=`echo $read | sed -e "s/_${readsuffix}\./_R${readsuffix}_001./g"`
            if [[ $verbose ]]; then
                echo " (file) $read contains $readkey"
            fi
        fi
        
        case $read in
            #check if contains lane before read
            *_L[0-9][0-9][0-9]_$readkey*)
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
                rename -f "s/_$readkey/_L001_$readkey/" ${read}*
                #update file variable
                read=`echo $read | sed -e "s/_${readkey}/_L001_${readkey}/g"`
                list[$j]=$read
            ;;
        esac
        case $read in
            #check if contains sample before lane
            *_S[0123456789]_L[0-9][0-9][0-9]*)
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
                k=$((${j} + 1))
                rename -f "s/_L([0123456789][0123456789][0123456789])/_S${k}_L\1/" ${read}*
                #update file variable
                read=`echo $read | sed -e "s/_L\([0123456789][0123456789][0123456789]\)/_S${k}_L\1/g"`
                list[$j]=$read
            ;;
        esac
            #check if contains sample before lane
            if [[ $read == *"_${readkey}_001."* ]]  || [[ $read == *"_${readkey}_001" ]]; then
                if [[ $verbose ]]; then
                    echo "  $read compatible with suffix"
                fi
            else
                #rename file
                if [[ $verbose ]]; then
                    echo "***Warning: file $read does not have suffix in its name. Suffix 001 is given.***"
                    echo "  renaming $read ..."
                fi
                if [[ -f $(find $(dirname ${read}) -name $(basename ${read})'*.gz') ]]; then
                    rename -f "s/_${readkey}.(.*).gz/_${readkey}_001\.\$1.gz/g" ${read}*gz
                fi
                if [[ ${read} == *.gz ]]; then
                    rename -f "s/_${readkey}.(.*).gz/_${readkey}_001\.\$1.gz/g" ${read}
                fi
                if [[ -f $(find $(dirname ${read}) -name $(basename ${read})'*.fastq') ]]; then
                    rename -f "s/_${readkey}\.(.*)/_${readkey}_001\.\$1/" ${read}*.fastq
                fi
                if [[ -f $(find $(dirname ${read}) -name $(basename ${read})'*.fq') ]]; then
                    rename -f "s/_${readkey}\.(.*)/_${readkey}_001\.\$1/" ${read}*.fq
                fi
                if [[ ${read} == *.fastq ]]; then
                    rename -f "s/_${readkey}*\.fastq/_${readkey}_001\.fastq/" ${read} ${read}.gz
               fi
              if [[ ${read} == *.fq ]]; then
                   rename -f "s/_${readkey}*\.fq/_${readkey}_001\.fq/" ${read}
              fi
              #update file variable
              if [[ ${read} == *.gz ]] || [[ ${read} == *.fastq ]] || [[ ${read} == *.fq ]] || [[ -f ${read} ]]; then
                  #assumes read name already contains . in file extension
                  read=`echo $read | sed -e "s/_${readkey}*\./_${readkey}_001\./g"`
              else
                  #replace everything after read key (R1, R2, I1, I2) with 001 suffix (detects file later)
                  rename -f "s/_${readkey}.*/_${readkey}_001/g" ${read}
                  read=`echo $read | sed -e "s/_${readkey}*/_${readkey}_001/g"`
              fi
              #remove characters after read key (R1, R2, I1, I2) required as above
              if [[ ${read} != *_${readkey}_001.* ]] && [[ ${read} != *.* ]]; then
                  rename -f "s/_${readkey}_*\./_${readkey}_001\./" ${read}
                  read=`echo $read | sed -e "s/_${readkey}*\./_${readkey}_001\./g"`
              elif [[ ${read} != *_${readkey}_001 ]] || [[ ${read} != *_${readkey}*00#1 ]]; then
                  rename -f "s/_${readkey}*_001/_${readkey}_001/" ${read#}
                  read=`echo $read | sed -e "s/_${readkey}*_001/_${readkey}_001/g"`
              fi
              if [[ ${read} == *_${readkey}_*_001.* ]]; then
                   rename -f "s/_${readkey}*_001\./_${readkey}_001\./" ${read}
                   read=`echo $read | sed -e "s/_${readkey}*_001\./_${readkey}_001\./g"`
              fi
              
              list[$j]=$read
        fi
        
        #allow detection of file extension (needed for --file input)
        if [[ -f ${read} ]] || [[ -h ${read} ]]; then
            echo "    $read file found"
        elif [[ -f ${read}.fq ]] || [[ -h ${read}.fq ]]; then
            read=${read}.fq
        elif [[ -f ${read}.fastq ]] || [[ -h ${read}.fastq ]]; then
            read=${read}.fastq
        elif [[ -f ${read}.fq.gz ]] || [[ -h ${read}.fq.gz ]]; then
            gunzip -f -k ${read}.fq.gz
            read=${read}.fq
        elif [[ -f ${read}.fastq.gz ]] || [[ -h ${read}.fastq.gz ]]; then
            gunzip -f -k ${read}.fastq.gz
            read=${read}.fastq.gz
        else
            echo "Error: $read not found"
            exit 1
        fi
            list[$j]=$read
    done
    
    if [[ $readkey == "R1" ]]; then
        read1=("${list[@]}")
    elif [[ $readkey == "R2" ]]; then
        read2=("${list[@]}")
    elif [[ $readkey == "R3" ]]; then
        read3=("${list[@]}")
    elif [[ $readkey == "R4" ]]; then
        read4=("${list[@]}")
    elif [[ $readkey == "I1" ]]; then
        index1=("${list[@]}")
    elif [[ $readkey == "I2" ]]; then
        index2=("${list[@]}")
    fi
done

if [[ $verbose ]]; then
    echo " Input files post-curation 3"
    echo "  ${#read1[@]}files - read1s: ${read1[@]}"
    echo "  ${#read2[@]}files - read2s: ${read2[@]}"
    if [[ ${#read3[@]} -gt 0 ]]; then
        echo "  ${#read3[@]}files - R3s: ${read3[@]}"
    fi
    if [[ ${#read4[@]} -gt 0 ]]; then
        echo "  ${#read4[@]}files - R4s: ${read4[@]}"
    fi
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]}files - I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]}files - I2s: ${index2[@]}"
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
        echo " checking for index I1 in R2 files ..."
        for ii in $(seq 1 1 ${#read1[@]}); do
             #iterate over read1 inputs
             indexfile=${read1[$(( $ii - 1 ))]}
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
        echo " checking for index I2 in R3 files ..."
        for ii in $(seq 1 1 ${#read1[@]}); do  
             #iterate over read1 inputs
             indexfile=${read1[$(( $ii - 1 ))]}
             #derive I2 filename from R1 filename
             indexfile=$(echo $indexfile | perl -pne 's/(.*)_R1/$1_R3/' )
             #only add index2 files to list variable if file exists (this will only run if I2 not found above)
             if [[ -f $indexfile ]] || [[ -f ${indexfile}.gz ]] || [[ -f $indexfile.fastq ]] || [[ -f ${indexfile}.fastq.gz ]] || [[ -f $indexfile.fq ]] || [[ -f ${indexfile}.fq.gz ]]; then
                 index2+=("$indexfile")
             fi
        done
    fi
    #checking for R1 and R4 index files
    echo " checking for read R4 files ..."
    for ii in $(seq 1 1 ${#read1[@]}); do
         #iterate over read1 inputs
         indexfile=${read1[$(( $ii - 1 ))]}
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
            echo "***WARNING: note that ${technology} expects dual indexes: I1 and I2 OR R2 and R3***"
        fi
    fi
fi


#generate missing indexes if required (generating I1 and I2)
if [[ "$technology" == "indrop-v3" ]] ||  [[ "$technology" == "icell8-full-length" ]] || [[ "$technology" == "sciseq2" ]] || [[ "$technology" == "sciseq3" ]] || [[ "$technology" == "scifiseq" ]] || [[ "$technology" == "smartseq2" ]] ||[[ "$technology" == "smartseq3" ]] || [[ "$technology" == "strt-seq-2i" ]] || [[ "$technology" == "strt-seq-2018" ]] || [[ "$technology" == "bravo" ]] || [[ "$technology" == "vasa-drop" ]]; then
   echo "dual indexes I1 and I2 required for $technology"
   if [[ ${#index2[@]} -le 0 ]]; then
       echo " automatically generating I1 and I2 index files from file headers"
       index1=("${read1[@]}")
       index2=("${read1[@]}")
       for ii in ${!read1[@]}; do
           #iterate over read1 inputs
           R1_file=${read1[$(( $ii - 1 ))]}
           R2_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_R2/' )
           I1_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_I1/' )
           I2_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_I2/' )
           
           if [[ $verbose ]]; then
               echo $R1_file
               echo $R2_file
               echo $I1_file
               echo $I2_file
           fi
           #copies index 1 to next line (1st to 2nd) and deletes 3rd line
           cat $R1_file | sed -E "s/ (.):(.):(.):(.*)\+(.*)$/ \1:\2:\3:\4+\5\n\4/g" | perl -n -e "print unless ($. % 5 == 3)" > $I1_file
           indexlength=$(($(head $I1_file -n 2 | tail -n 1 | wc -c) -1))
           qualscores=$(seq 1 $indexlength | xargs -I {} printf I)
           if [[ $verbose ]]; then
               echo index of length $indexlength gives quality score $qualscores
           fi
           perl -pni -e "s/^.*$/${qualscores}/g if ($. % 4 == 0)" $I1_file
           #copies index 2 to next line (1st to 2nd) and deletes 3rd line
           cat $R1_file | sed -E "s/ (.):(.):(.):(.*)\+(.*)$/ \1:\2:\3:\4+\5\n\5/g" | perl -n -e "print unless ($. % 5 == 3)" >  $I2_file
           index2length=$(($(head $I2_file -n 2 | tail -n 1 | wc -c) - 1))
           qualscores2=$(seq 1 $index2length | xargs -I {} printf I)
           if [[ $verbose ]]; then
               echo index2 of length $index2length gives quality score $qualscores2
           fi
           perl -pni -e "s/^.*$/${qualscores2}/g if ($. % 4 == 0)" $I2_file
           index1+=("$I1_file")
           index2+=("$I2_file")
        done
        if [[ $verbose ]]; then
            echo index1: $index1
            echo index2: $index2
        fi
    else
        echo " dual indexes found"
    fi
fi

if [[ "$technology" == "quartz-seq" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "strt-seq-c1" ]]; then
    echo "dual indexes I1 and I2 required for $technology"
    if [[ ${#index2[@]} -le 1 ]]; then
        echo " automatically generating I1 index files from file headers"
        index1=("${read1[@]}")
        for ii in ${!read1[@]}; do
            #iterate over read1 inputs
            R1_file=${read1[$(( $ii - 1 ))]}
            R2_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_R2/' )
            I1_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_I1/' )
            if [[ $verbose ]]; then
                echo $R1_file
                echo $R2_file
                echo $I1_file
            fi
            #copies index 1 to next line (1st to 2nd) and deletes 3rd line (only if index 1 doesn't contain '+' character)
            cat $R1_file | sed -E "/x/! s/ (.):(.):(.):(.*)$/ \1:\2:\3:\4$\n\4/g" > $I1_file
            linediff=$(grep -n "^+" $I1_file | head -n 2 | cut -d":" -f 1 | awk 'NR==1{p=$1;next} END{print $1-p}')
            if [[ $linediff -eq 5 ]];then
                #remove lines if matched only
                perl -n -e "print unless ($. % 5 == 3)" > $I1_file
            else
                cat $R1_file | sed -E "s/ (.):(.):(.):(.*)\+(.*)$/ \1:\2:\3:\4+\5\n\4/g" | perl -n -e "print unless ($. % 5 == 3)" > $I1_file
            fi
            indexlength=$(($(head $I1_file -n 2 | tail -n 1 | wc -c) - 1))
            qualscores=$(seq 1 $indexlength | xargs -I {} printf I)
            if [[ $verbose ]]; then
                echo index of length $indexlength gives quality score $qualscores
            fi
            perl -pni -e "s/^.*$/${qualscores}/g if ($. % 4 == 0)" $I1_file
            index1+=("$I1_file")
        done
        if [[ $verbose ]]; then
            echo index1: $index1
        fi
    else
        echo " index found"
    fi
fi

#inverting R1 and R2 for specific technologies
if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "splitseq" ]] || [[ "$technology" == "splitseq2" ]] || [[ "$technology" == "strt-seq-2018" ]] || [[ "$technology" == "vasa-drop" ]]; then
    #invert read1 and read2
    echo "***WARNING: technology is set to ${technology}. barcodes on Read 2 will be used***"
    tmp=$read1
    read1=$read2
    read2=$tmp
    tmp=""
fi

if [[ $verbose ]]; then
    echo " Input files post-curation 4"
    echo "  ${#read1[@]}files - read1s: ${read1[@]}"
    echo "  ${#read2[@]}files - read2s: ${read2[@]}"
    if [[ ${#read3[@]} -gt 0 ]]; then
        echo "  ${#read3[@]}files - R3s: ${read3[@]}"
    fi
    if [[ ${#read4[@]} -gt 0 ]]; then
        echo "  ${#read4[@]}files - R4s: ${read4[@]}"
    fi
    if [[ ${#index1[@]} -gt 0 ]]; then
        echo "  ${#index1[@]}files - I1s: ${index1[@]}"
    fi
    if [[ ${#index2[@]} -gt 0 ]]; then
        echo "  ${#index2[@]}files - I2s: ${index2[@]}"
    fi
    echo "  input files adjusted for technology-specific conditions"
fi
##########



#####Input file curation 5: Catpuring sample name#####
#checking the quality of fastq file names
read12=("${read1[@]}" "${read2[@]}")
if [[ ${#read3[@]} -gt 0 ]] && [[ ${#read4[@]} -gt 0 ]]; then
    read12=("${read1[@]}" "${read2[@]}" "${read3[@]}" "${read4[@]}")
elif [[ ${#index1[@]} -gt 0 ]] && [[ ${#read3[@]} -gt 0 ]]; then
    read12=("${read1[@]}" "${read2[@]}" "${read3[@]}" "${index1[@]}")
elif [[ ${#index2[@]} -gt 0 ]]; then
     read12=("${read1[@]}" "${read2[@]}" "${index1[@]}" "${index2[@]}")
elif [[ ${#index1[@]} -gt 0 ]]; then
    read12=("${read1[@]}" "${read2[@]}" "${index1[@]}")
fi

if [[ $verbose ]]; then
    echo "files to be curated:" ${read12[@]}
fi

for fq in "${read12[@]}"; do
    if [[ $verbose ]]; then
        echo " read file: $fq"
    fi
    name=`basename $fq`
    name=${name%.*}
    fields=`echo ${name} | grep -o "_" | wc -l`
    fields=$(($fields + 1))
    if [[ $fields -le 4 ]]; then
        name_fields=$(($fields + 1))
    else
        name_fields=$fields
    fi
    sn=`echo ${name} | cut -f 1-$((${name_fields} - 4)) -d'_'`
    # removes leading zeroes from lane number
    lane=`echo ${name} | cut -f $((${fields} - 2)) -d'_' | sed -E 's/L([123456789][0123456789][0123456789])/\1/' | sed -E 's/L0([123456789][0123456789])/\1/' | sed -E 's/L00([0123456789])/\1/'`
    # sets lane to 0 if none found
    if [[ -z ${lane} ]]; then
        lane=0
    fi
    LANE+=($lane)
    if [[ ${fields} -le 4 ]]; then
        echo "***WARNING: filename $fq is not following the naming convention. (e.g. EXAMPLE_S1_L001_R1_001.fastq)***";
        #exit 1
    elif [[ $fq != *'.fastq'* ]] && [[ $fq != *'.fq'* ]]; then
        echo "Error: $fq does not have a .fq or .fastq extention."
        exit 1
    elif [[ ${sn} =~ "." ]]; then
        echo "Error: $fq has a period \".\" within its sample name. Remove it to run Cell Ranger."
        exit 1
    fi
    
    if [[ $verbose ]]; then
        echo "  $sn(extracted from file) <- $fq"
        echo "  SAMPLE NAME: $SAMPLE"
    fi
    
    if [[ ${sn} != $SAMPLE ]]; then
        if [[ -z $SAMPLE ]]; then
            if [[ $verbose ]]; then
                echo "  setting SAMPLE NAME to ${sn}"
            fi
            SAMPLE=${sn}
        else
            echo "Error: some samples are labeled $SAMPLE while others are labeled $sn. cellranger can only handle files from one sample at a time."
            exit 1
        fi
    fi
done

if [[ $verbose ]]; then
     echo "read1, read2, index1, and index2 file curation complete"
fi

LANE=$(echo "${LANE[@]}" | tr ' ' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
##########



#####Set the input barcode file#####
if [[ $verbose ]]; then
    echo " setting whitelist barcode file."
fi

custombarcodes=false
if [[ -n "$barcodefile" ]]; then
    if [[ ! -f $barcodefile ]]; then
        echo "Error: File selected for option --barcodefile does not exist"
        exit 1
    else
        #getting absolute path
        barcodefile=$(readlink -f $barcodefile)
        custombarcodes=true
        #allowing WellList from ICELL8 and other well-based techniques
        if [[ "$technology" == "bd-rhapsody" ]] || [[ "$technology" == "bd-rhapsody-v2" ]] || [[ "$technology" == "bravo" ]] || [[ "$technology" == "celseq" ]] || [[ "$technology" == "celseq2" ]] || [[ "$technology" == "fluidigm-c1" ]] || [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "icell8" ]] || [[ "$technology" == "quartz-seq" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]] || [[ "$technology" == "quartz-seq2*" ]] || [[ "$technology" == "microwellseq" ]] || [[ "$technology" == "smartseq*" ]] || [[ "$technology" == "seqwell" ]] || [[ "$technology" == "sciseq2" ]] || [[ "$technology" == "sciseq3" ]] || [[ "$technology" == "scifiseq" ]] || [[ "$technology" == "splitseq" ]] || [[ "$technology" == "splitseq2" ]] || [[ "$technology" == "pip-seq-v0" ]] || [[ "$technology" == "pip-seq-v1" ]] || [[ "$technology" == "pip-seq-v2" ]]  || [[ "$technology" == "pip-seq-v3" ]]  || [[ "$technology" == "pip-seq-v4" ]] || [[ "$technology" == "vasa-plate" ]] || [[ "$technology" == "custom" ]]; then
            seg=$'\t'
            n_col=$(awk -F'\t' '{print NF}' $barcodefile | sort -nu | tail -n 1)
            if [[ $n_col -eq 1 ]]; then
                seg=","
                n_col=$(awk -F',' '{print NF}' $barcodefile | sort -nu | tail -n 1)
            fi
            
            if [[ $n_col -gt 1 ]]; then
                new_barcodefile=${barcodefile%.*}_barcode.txt
                
                #get column with barcodes
		col_n=$(head -n 1 $barcodefile | tr "${seg}" "\n" | grep -n "[Bb]arcode" | head -n 1 | cut -d":" -f 1)
                if [[ -z $col_n ]]; then
                    col_n=1
                fi
                
                if [[ $(head -n 1 $barcodefile | grep "[Bb]arcode") ]]; then
                    #removes header (1st line) containing colname "barcode"
                    tail -n $(($(wc -l $barcodefile | cut -d" " -f 1) - 1)) $barcodefile | cut -f $col_n -d "$seg" > ${new_barcodefile}
                else
                    echo "***WARNING: barcode file has multiple columns and none named 'barcode(s)'***"
                    echo "Note: please check that column 1 of the barcode file contains the barcodes as required"
                    head $barcodefile | cut -f 1 -d "${seg}"
                    #assumes no headers
                    cut -f 1 -d "${seg}" $barcodefile > ${new_barcodefile}
                fi
                barcodefile=${new_barcodefile}
            fi
        fi
    fi
else
    if [[ "$technology" == "10x"* ]]; then
        barcodefile="default:10x"
    elif [[ "$technology" == "bd-rhapsody" ]] || [[ "$technology" == "bd-rhapsody-v2" ]]; then
        barcodefile=${whitelistdir}/bd_rhapsody_barcode.txt
        if [[ ! -f ${whitelistdir}/bd_rhapsody_barcode.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "bravo" ]]; then
        barcodefile=${whitelistdir}/KAPA_UDI_dual_barcodes.txt
        if [[ ! -f ${whitelistdir}/KAPA_UDI_dual_barcodes.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "celseq" ]]; then
        barcodefile=${whitelistdir}/bc_celseq1.txt
        if [[ ! -f ${whitelistdir}/bc_celseq1.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "celseq2" ]]; then
        barcodefile=${whitelistdir}/bc_celseq2.txt
        if [[ ! -f ${whitelistdir}/bc_celseq2.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "fluidigm-c1" ]] || [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]]; then
        barcodefile=${whitelistdir}/Illumina_Nextera_dual_barcodes.txt
        if [[ ! -f ${whitelistdir}/Illumina_Nextera_dual_barcodes.txt ]]; then
            echo "  generating combination of I1 and I2 barcodes ..."
        fi
    elif [[ "$technology" == "icell8" ]]; then
        barcodefile=${whitelistdir}/ICELL8_barcode.txt
        echo "***WARNING: selected barcode file (${barcodefile}) contains barcodes for all wells in ICELL8. valid barcode will be an overestimate***"
    elif [[ "$technology" == "icell8-5-prime" ]]; then
        barcodefile=${whitelistdir}/ICELL8_TCR_barcode.txt
    elif [[ "$technology" == "icell8-full-length" ]]; then
        barcodefile=${whitelistdir}/SmartSeq_ICELL8_dual_barcodes.txt
        if [[ ! -f ${whitelistdir}/SmartSeq_ICELL8_dual_barcodes.txt ]]; then
            echo "  generating combination of I1 and I2 barcodes ..."
        fi
    elif [[ "$technology" == "marsseq-v2" ]]; then
        barcodefile=${whitelistdir}/MARS-Seq2_barcode.txt
    elif [[ "$technology" == "microwellseq" ]]; then
        barcodefile=${whitelistdir}/microwellseq_barcode.txt
        if [[ ! -f ${whitelistdir}/microwellseq_barcode.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "quartz-seq" ]]; then
        indexlength=$(($(head $index1([0]) -n 2 | tail -n 1 | wc -c) -1))
        if [[ -f $(echo ${index2}([0])) ]]; then
            index2length=$(($(head $index2([0]) -n 2 | tail -n 1 | wc -c) -1))
            barcodelength=$(($indexlength+$index2length))
            if [[ -f ${whitelistdir}/Illumina_dual_barcodes.txt ]]; then
                cat ${whitelistdir}/Illumina_TruSeq_Index1_i7_barcodes.txt ${whitelistdir}/Illumina_Nextera_Index1_i7_barcodes.txt | sort | uniq > ${whitelistdir}/Illumina_Index1_i7_barcodes.txt
                cat ${whitelistdir}/Illumina_TruSeq_Index2_i5_barcodes.txt ${whitelistdir}/Illumina_Nextera_Index2_i5_barcodes.txt | sort | uniq > ${whitelistdir}/Illumina_Index2_i5_barcodes.txt
                join -j 9999 ${whitelistdir}/Illumina_Index1_i7_barcodes.txt ${whitelistdir}/Illumina_Index2_i5_barcodes.txt | sed "s/ //g" > ${whitelistdir}/Illumina_dual_barcodes.txt
            fi
            barcodefile=${whitelistdir}/Illumina_dual_barcodes.txt
        else
            barcodelength=$indexlength
            if [[ $indexlength -eq 6 ]]; then
                barcodefile=${whitelistdir}/Illumina_TruSeq_LT_Index1_i7_barcodes.txt
            else
                cat ${whitelistdir}/Illumina_TruSeq_Index1_i7_barcodes.txt ${whitelistdir}/Illumina_Nextera_Index1_i7_barcodes.txt >${whitelistdir}/Illumina_Index1_i7_barcodes.txt
                barcodefile=${whitelistdir}/Illumina_Nextera_Index1_i7_barcodes.txt
            fi
        fi
    elif [[ "$technology" == "quartz-seq2-384" ]]; then
        barcodefile=${whitelistdir}/Quartz-Seq2-384_barcode.txt
    elif [[ "$technology" == "quartz-seq2-1536" ]]; then
        barcodefile=${whitelistdir}/Quartz-Seq2-1536_barcode.txt
    elif [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]]; then
        echo "***WARNING: whitelist for ${technology} is modified from the original barcodes (https://github.com/indrops/indrops/tree/master/ref/barcode_lists), first 8 bp of list1 and list2 are joind to generate a 16 bp barcode***"
        barcodelength=${minlength}
        if [[ $verbose ]]; then
            echo "  barcode adjusted to ${barcodelength} bp to match the length in the default whitelist for ${technology}"
        fi
        if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
            barcodefile=${whitelistdir}/inDrop-v1_barcodes.txt
        elif [[ "$technology" == "indrop-v3" ]]; then
            barcodefile=${whitelistdir}/inDrop-v3_barcodes.txt
            echo "***WARNING: combination of list1 and list2 from indrop-v2 (https://github.com/indrops/indrops/issues/32)***"  
        fi
    elif [[ "$technology" == "pip-seq-v0" ]]; then
        barcodefile=${whitelistdir}/pip-seq-v0_barcodes.txt
        if [[ ! -f ${whitelistdir}/pip-seq-v0_barcodes.txt ]]; then
            echo "  generating combination of barcodes ..."
        fi
    elif [[ "$technology" == "pip-seq-v1" ]]; then
        barcodefile=${whitelistdir}/pip-seq-v1_barcodes.txt
        if [[ ! -f ${whitelistdir}/pip-seq-v1_barcodes.txt ]]; then
            echo "  generating combination of barcodes ..."
        fi
    elif [[ "$technology" == "pip-seq-v2" ]]; then
        barcodefile=${whitelistdir}/pip-seq-v2_barcodes.txt
        if [[ ! -f ${whitelistdir}/pip-seq-v2_barcodes.txt ]]; then
            echo "  generating combination of barcodes ..."
        fi
    elif [[ "$technology" == "pip-seq-v3" ]] || [[ "$technology" == "pip-seq-v4" ]]; then
        barcodefile=${whitelistdir}/pip-seq-v3_barcodes.txt
        if [[ ! -f ${whitelistdir}/pip-seq-v3_barcodes.txt ]]; then
            echo "  generating combination of barcodes ..."
        fi
    elif [[ "$technology" == "sciseq2" ]]; then
            barcodefile=${whitelistdir}/sciseq2_barcode.txt
            if [[ ! -f ${whitelistdir}/sciseq2_barcode.txt ]]; then
                echo "  generating combination of I1, I2, and RT barcodes ..."
            fi
    elif [[ "$technology" == "sciseq3" ]]; then
            barcodefile=${whitelistdir}/sciseq3_barcode.txt
            if [[ ! -f ${whitelistdir}/sciseq3_barcode.txt ]]; then
                echo "  generating combination of I1, I2, and RT barcodes ..."
            fi
    elif [[ "$technology" == "scifiseq" ]]; then
            barcodefile=${whitelistdir}/scifi-seq_barcode.txt
            if [[ ! -f ${whitelistdir}/scifi-seq_barcode.txt ]]; then
                echo "  generating combination of I1, I2, and RT barcodes ..."
            fi
    elif [[ "$technology" == "splitseq" ]] || [[ "$technology" == "splitseq2" ]]; then
            barcodefile=${whitelistdir}/splitseq_3_barcodes.txt
            if [[ ! -f ${whitelistdir}/splitseq_3_barcodes.txt ]]; then
                echo "  generating combination of I1, I2, and RT barcodes ..."
            fi
    elif [[ "$technology" == "smartseq2" ]]; then
        barcodefile=${whitelistdir}/SmartSeq2_full_barcodes.txt
    elif [[ "$technology" == "smartseq3" ]]; then
        barcodefile=${whitelistdir}/SmartSeq3_full_barcodes.txt
    elif [[ "$technology" == "strt-seq" ]]; then
        barcodefile=${whitelistdir}/STRTSeq_barcode.txt
    elif [[ "$technology" == "strt-seq-c1" ]]; then
        barcodefile=${whitelistdir}/STRTSeqC1_barcode.txt
    elif [[ "$technology" == "strt-seq-2i" ]]; then
        barcodefile=${whitelistdir}/STRTSeq2i_barcode.txt
    elif  [[ "$technology" == "strt-seq-2018" ]]; then
        barcodefile=${whitelistdir}/strt_fan2018_barcode_96_8bp.txt
    elif [[ "$technology" == "vasa-plate" ]]; then
        barcodefile=${whitelistdir}/bc_vasa_plate.txt
        if [[ ! -f ${whitelistdir}/bc_vasa_plate.txt ]]; then
            echo "  generating combination of I1, I2, and RT barcodes ..."
        fi
    elif [[ "$technology" == "vasa-drop" ]]; then
        if [[ $custombarcodes == false ]]; then
            echo "WARNING: For VASA-drop, you are recommended produce the barcode file for each library with --barcodefile as input parameter."
            barcodefile=${whitelistdir}/bc_vasa_drop.txt
            if [[ ! -f ${whitelistdir}/bc_vasa_drop.txt ]]; then
                echo "  generating combination of I1, I2, and RT barcodes ..."
            fi
        else
            if [[ ! -f $barcodefile ]]; then
                echo "Error: File selected for option --barcodefile does not exist"
                exit 1
            fi
        fi
    else
        echo "***WARNING: whitelist for ${technology} will be all possible combinations of ${minlength} bp. valid barcode will be 100% as a result***"
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
        echo "  barcodefile file exists: ${barcodefile}"
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
        if [[ "$technology" == "bd-rhapsody" ]] || [[ "$technology" == "bd-rhapsody-v2" ]]; then
            if [[ ! -f ${whitelistdir}/bd_rhapsody_barcode.txt ]]; then
                #generates all combinations of I1-I2-R1 barcodes
                join -j 9999 ${whitelistdir}/bd_rhapsody_cell_label_section1.txt ${whitelistdir}/bd_rhapsody_cell_label_section2.txt | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/bd_rhapsody_cell_label_section3.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/bd_rhapsody_barcode.txt
            fi
        elif [[ "$technology" == "bravo" ]]; then
            if [[ ! -f barcodefile=${whitelistdir}/KAPA_UDI_dual_barcodes.txt ]]; then
                #generates all combinations of I1-I2-R1 barcodes
                join -j 9999 ${whitelistdir}/KAPA_UDI_Index1_i7.txt ${whitelistdir}/KAPA_UDI_Index5_i5.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/KAPA_UDI_dual_barcodes.txt
            fi
        elif [[ "$technology" == "fluidigm-c1" ]] || [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]] || [[ "$technology" == "smartseq2" ]] || [[ "$technology" == "smartseq3" ]]; then
            if [[ ! -f ${whitelistdir}/Illumina_Nextera_dual_barcodes.txt ]];then
                #generates all combinations of I1-I2 barcodes
                join -j 9999 ${whitelistdir}/Illumina_Nextera_Index1_i7_barcodes.txt ${whitelistdir}/Illumina_Nextera_Index2_i5_barcodes.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/Illumina_Nextera_dual_barcodes.txt
            fi
        elif [[ "$technology" == "icell8-full-length" ]]; then
            if [[ ! -f ${whitelistdir}/SmartSeq_ICELL8_dual_barcodes.txt ]]; then
                #generates all combinations of I1-I2 barcodes
                join -j 9999 ${whitelistdir}/ICELL8_full_length_Index1_i7_barcodes.txt ${whitelistdir}/ICELL8_full_length_Index2_i5_barcodes.txt | sed "s/ //g" \
                sort | uniq \
                > ${whitelistdir}/SmartSeq_ICELL8_dual_barcodes.txt
            fi
        elif [[ "$technology" == "indrop-v"* ]]; then
            if [[ "$technology" == "indrop-v1" ]] || [[ $technology"" == "indrop-v2" ]]; then
                sed -E 's/.*(.{8})/\1/g' ${whitelistdir}/inDrop_gel_barcode1_list_revcomp.txt > ${whitelistdir}/inDrop_gel_barcode1_list_revcomp_tail.txt
                perl ${MAKEINDROPBARCODES} ${whitelistdir}/inDrop_gel_barcode1_list_revcomp_tail.txt ${whitelistdir}/inDrop_gel_barcode2_list_revcomp.txt v1 ${whitelistdir}
            elif [[ "$technology" == "indrop-v3" ]]; then
                #allow for barcodes in index (I1) and R1
                perl ${MAKEINDROPBARCODES} ${whitelistdir}/inDrop_gel_barcode1_list.txt ${whitelistdir}/inDrop_gel_barcode2_list.txt v3 ${whitelistdir}
            fi
        elif [[ "$technology" == "microwellseq" ]]; then
            if [[ ! -f ${whitelistdir}/microwellseq_barcode.txt ]]; then
                #generates all combinations of R1 barcodes
                join -j 9999 ${whitelistdir}/microwell-seq_barcodeA.txt ${whitelistdir}/microwell-seq_barcodeB.txt | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/microwell-seq_barcodeC.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/microwellseq_barcode.txt
            fi
        elif [[ "$technology" == "pip-seq-v0" ]]; then
            if [[ ! -f ${whitelistdir}/pip-seq-v0_barcodes.txt ]]; then
                #generates all combinations of R1 barcodes
                join -j 9999 ${whitelistdir}/pip-seq_v0_bc1.tsv ${whitelistdir}/pip-seq_v0_bc2.tsv | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/pip-seq_v0_bc3.tsv | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/pip-seq-v0_barcodes.txt
            fi
        elif [[ "$technology" == "pip-seq-v1" ]]; then
            if [[ ! -f ${whitelistdir}/pip-seq-v1_barcodes.txt ]]; then
                #generates all combinations of R1 barcodes
                join -j 9999 ${whitelistdir}/pip-seq_v1_bc1.tsv ${whitelistdir}/pip-seq_v1_bc1.tsv | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/pip-seq-v1_barcodes.txt
            fi
        elif [[ "$technology" == "pip-seq-v2" ]]; then
            if [[ ! -f ${whitelistdir}/pip-seq-v2_barcodes.txt ]]; then
                #generates all combinations of R1 barcodes
                join -j 9999 ${whitelistdir}/pip-seq_v2_bc1.tsv ${whitelistdir}/pip-seq_v2_bc2.tsv | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/pip-seq_v2_bc3.tsv | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/pip-seq-v2_barcodes.txt
            fi
        elif [[ "$technology" == "pip-seq-v3" ]] | [[ "$technology" == "pip-seq-v4" ]]; then
            if [[ ! -f ${whitelistdir}/pip-seq-v4_barcodes.txt ]]; then
                #generates all combinations of R1 barcodes
                join -j 9999 ${whitelistdir}/pip-seq_v3_bc1.tsv ${whitelistdir}/pip-seq_v2_bc3.tsv | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/pip-seq_v3_bc3.tsv | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/pip-seq_v3_bc4.tsv | sed "s/ //g" |
                sort | uniq \
                > ${whitelistdir}/pip-seq-v3_barcodes.txt
            fi
        elif [[ "$technology" == "sciseq2" ]]; then
            if [[ ! -f ${whitelistdir}/sciseq2_barcode.txt ]]; then
                #generates all combinations of I1-I2-R2 barcodes
                join -j 9999 ${whitelistdir}/sci-seq3_i7_barcodes.txt ${whitelistdir}/sci-seq3_i5_barcodes.txt | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/sci-seq3_rt_barcodes.txt | sed "s/ //g" | awk '!a[$0]++' | \
                sort | uniq \
                > ${whitelistdir}/sciseq2_barcode.txt
            fi
        elif [[ "$technology" == "sciseq3" ]]; then
            if [[ ! -f ${whitelistdir}/sciseq3_barcode.txt ]]; then
                #generates all combinations of I1-I2-R1 barcodes
                join -j 9999 ${whitelistdir}/sci-seq3_i7_barcodes.txt ${whitelistdir}/sci-seq3_i5_barcodes.txt | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/sci-seq3_hp_barcodes.txt | sed "s/ //g" | join -j 9999 - ${whitelistdir}/sci-seq3_rt_barcodes.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/sciseq3_barcode.txt
                ## to filter unique lines: awk '!a[$0]++' > ${whitelistdir}/sciseq3_barcode.txt
            fi
        elif [[ "$technology" == "scifiseq" ]]; then
            if [[ ! -f ${whitelistdir}/scifi-seq_barcode.txt ]]; then
                #generates all combinations of I1-I2-R1 barcodes
                join -j 9999 ${whitelistdir}/10x_atac_barcodes.txt ${whitelistdir}/ scifi-seq_rt_barcode.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/scifi-seq_barcode.txt
            fi
        elif [[ "$technology" == "splitseq" ]]; then
            if [[ ! -f ${whitelistdir}/splitseq_3_barcodes.txt ]]; then
                #generates all combinations of I1-I2-R2 barcodes
                join -j 9999 ${whitelistdir}/split-seq_round1_barcodes.txt ${whitelistdir}/split-seq_round2_barcodes.txt | sed "s/ //g" | \
                join -j 9999 - ${whitelistdir}/split-seq_round3_barcodes.txt | sed "s/ //g" | awk '!a[$0]++' | \
                sort | uniq \
                > ${whitelistdir}/splitseq_3_barcodes.txt
            fi
        elif [[ "$technology" == "strt-seq-2i" ]]; then
            if [[ ! -f ${whitelistdir}/STRTSeq2i_barcode.txt ]]; then
                #generates all combinations of I1-I2-R1 barcodes
                join -j 9999 ${whitelistdir}/AllPossibilities_5_barcodes.txt ${whitelistdir}/STRTSeqC1_barcode.txt | sed "s/ //g" | \
                sort | uniq \
                > ${whitelistdir}/STRTSeq2i_barcode_barcode.txt
           fi
        else
            #generating permutations of ATCG of barcode length (non-standard evaluation required to run in script)
            if [[ ${barcodelength} -ge 12 ]]; then
                echo "  ... generating all permutations of A,T,C,G of length ${barcodelength}"
                echo "***WARNING: for large barcodes this could take a lot of time and memory***"
                echo "  Please use a known barcode whitelist if possible"
            fi
            echo $(eval echo $(for ii in $(eval echo {1..${barcodelength}}); do echo "{A,T,C,G}"; done | tr "\n" " " | sed "s/ //g" | xargs -I {} echo {})) | sed 's/ /\n/g' | sort | uniq > ${barcodefile}
        fi
    fi
fi

if [[ $verbose ]]; then
    echo " barcodefile generated"
fi
##########



####check if reference is present#####
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
    if [[ $umilength -ge 11 ]]; then
        if [[ "$chemistry" == "SC3Pv1" ]] || [[ "$chemistry" == "SC3Pv2" ]]; then
            echo "Using 10x version 3 chemistry to support longer UMIs"
            chemistry="SC3Pv3"
            umi_default=12
        fi
    elif [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R1" ]] && [[ "$chemistry" != "SC5P-R2" ]] && [[ "$chemistry" != "fiveprime" ]]; then
        #use SC3Pv2 (umi length 10)
        echo "Using 10x version 2 chemistry to support UMIs"
        chemistry="SC3Pv2"
        umi_default=10
    fi
    if [[ "$chemistry" != "threeprime" ]] && [[ "$chemistry" != "fiveprime" ]] && [[ "$chemistry" != "SC3Pv1" ]] && [[ "$chemistry" != "SC3Pv2" ]] && [[ "$chemistry" != "SC3Pv3" ]] && [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R1" ]] && [[ "$chemistry" != "SC5P-R2" ]]; then
       echo "Error: option --chemistry must be SC3Pv3, SC3Pv2, SC5P-PE, SC5P-R1, or SC5P-R2"
       exit 1
    fi
fi
if [[ "$technology" == "10x" ]]; then
    #use SC3Pv3 (umi length 12)
    if [[ "$chemistry" == "SC3Pv1" ]]; then
        echo "Accepted chemistry: $chemistry"
        barcode_default=14
        umi_default=10
    elif [[ "$chemistry" == "SC3Pv2" ]]; then
        echo "Accepted chemistry: $chemistry"
        barcode_default=16
        umi_default=10
    elif [[ "$chemistry" == "SC3Pv3" ]]; then
        echo "Accepted chemistry: $chemistry"
        barcode_default=16
        umi_default=12
    elif [[ "$chemistry" != "auto" ]]; then
        #use automatic chemistry detection
        echo "Detecting 10x chemistry automatically"
        chemistry="auto"
        #do not convert UMI
        umi_default=12
    fi
    # disable conversion for 10x
    barcodelength=${barcode_default}
    barrcodeadjust=0
    umilength=${umi_default}
    umiadjust=0
fi
if [[ "$technology" == "smartseq3" ]] || [[ "$technology" == "icell8-5-prime" ]]; then
    if [[ $verbose ]]; then
        echo "  Using $chemistry for $technology"
    fi
    if [[ "$chemistry" == "fiveprime" ]]; then
       chemistry="SC5P-PE"
    fi
    if [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R1" ]] && [[ "$chemistry" != "SC5P-R2" ]]; then
        if [[ $nonUMI == "false" ]]; then
            echo "Error: option --chemistry must be SC5P-PE, SC5P-R1 or SC5P-R2 for ${technology} (umi-based)"
            exit 1
        fi
    fi
fi
if [[ "$technology" == "smartseq2"  || ( "$technology" == "smartseq3" && $nonUMI == "true" ) ]]; then
    if [[ $verbose ]]; then
        echo "  Using $chemistry for $technology"
    fi
    if [[ "$chemistry" == "fiveprime" ]]; then
       chemistry="SC5P-PE"
    fi
    if [[ "$chemistry" == "SC5P-R1" ]]; then
        echo "  accurately mapping 5\' ends (filtering by tag sequence)"
    else
        echo "  mapping all reads (tag sequence removed)"
    fi
    if [[ "$chemistry" == "SC5P-PE" ]]; then
        echo "  mapping paired ends..."
    else
        echo "  mapping single-ends..."
    fi
    if [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R1" ]] && [[ "$chemistry" != "SC5P-R2" ]]; then
         echo "  using full-length sequences for read counts (default for ${technology})"
    fi
fi
##########



#####checking if jobmode matches expected input#####
if [[ "$jobmode" != "local" ]] && [[ "$jobmode" != "sge" ]] && [[ "$jobmode" != "lsf" ]] && [[ "$jobmode" != "slurm" ]] && [[ "$jobmode" != *"template" ]]; then
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
    echo "***WARNING: conversion was turned on because directory $crIN was not found***"
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
            currentbarcode=$(echo $barcodefile | rev | cut -d"/" -f 1 | rev)
            lastbarcode=$(echo $lastcall_p | rev | cut -d"/" -f 1 | rev)
            if [[ $verbose == true ]]; then
                echo "current file: $currentbarcode"
                echo "last file: $lastbarcode"
            fi
            if [[ ${barcodelength} == ${lastcall_b} ]] && [[ ${umilength} == ${lastcall_u} ]] && [[ ${barcodefile} == ${lastcall_p} ]]; then
                echo " call accepted: no conflict detected with other jobs currently running"
                #add current job to lock
                lock=$(($lock + 1))
                if [[ $setup == "false" ]]; then 
                    echo $lock > $lockfile
                fi
            #compare filenames (not paths) if not custom barcodes
            elif [[ ${barcodelength} == ${lastcall_b} ]] && [[ ${umilength} == ${lastcall_u} ]] && [[ $currentbarcode == "AllPossibilities_"* ]] && [[ $lastbarcode == "AllPossibilities_"* ]]; then
                echo " call accepted: no conflict detected with other jobs currently running"
                #add current job to lock
                lock=$(($lock + 1))
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
    echo "***WARNING: launch_universc.sh will exit once whitelist is converted***"
fi
echo "FORMAT: $technology"
if [[ $technology == "nadia" ]]; then
    echo "***WARNING: whitelist is converted for compatibility with $technology, valid barcodes cannot be detected accurately with this technology***"
fi
echo "BARCODES: ${barcodefile}"
if [[ ${#read1[@]} -eq 0 ]] && [[ ${#read1[@]} -eq 0 ]]; then
    echo "***WARNING: no FASTQ files were selected, launch_universc.sh will exit after setting up the whitelist***"
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
    echo "***WARNING: no description given, setting to ID value***"
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
    echo "***WARNING: --jobmode \"sge\" is recommended if running script with qsub***"
fi
echo "CONVERSION: $convert"
if [[ $convert == "false" ]]; then
    echo "***WARNING: adjustment for barcode and UMI length was skipped***"
fi
echo "##########"
echo ""
##########



####setup whitelist#####
if [[ $verbose ]]; then
    echo "lock $lock"
fi
if [[ $lock -eq 0 ]]; then
    echo "whitelist setup begin"
    echo "updating barcodes in $barcodedir for Cell Ranger version ${cellrangerversion} installed in ${cellrangerpath} ..."
    
    cd $barcodedir
    
    #restore assert functions if cellranger version is 3 or greater
    echo " restoring Cell Ranger"
    if [[ $verbose ]]; then
        echo "last call: $lastcall_p"
    fi
    if [[ $(echo "${cellrangerversion} 3.0.0" | tr " " "\n" | sort -V | tail -n 1) == ${cellrangerversion} ]]; then
        if [[ $verbose ]]; then
            echo "cellranger version 3.0.0 or greater"
        fi
        if [[ ! -z $barcodefile ]] && [[ $verbose ]]; then
            echo "barcodefile: $barcodefile"
        fi
        if [[ $technology == "10x" ]] && [[ $barcodefile == "default:10x" ]]; then
            #restore checking barcodes
            if [[ $verbose ]]; then
                echo "restore barcode checks"
            fi
            sed -i "s/#*#if gem_group == prev_gem_group/if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#*#assert barcode_idx >= prev_barcode_idx/assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/#*#assert np.array_equal(in_mc.get_barcodes(), barcodes)/assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        elif [[ $lastcall_p == "default:10x" ]] || [[ ! -f $lastcallfile ]]; then
            #disable checking barcodes
            if [[ $verbose ]]; then
                 echo "disable barcode checks"
            fi
            sed -i "s/if gem_group == prev_gem_group/#if gem_group == prev_gem_group/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${cellrangerpath}-cs/${cellrangerversion}/mro/stages/counter/report_molecules/__init__.py
            sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/molecule_counter.py
        fi
    fi
    
    #backup 10x navbar
    if [[ ! -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html ]];then
        cp ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.html ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html
        cp ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.html ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.temp.html
    fi
    if [[ $technology == "10x" ]] && [[ $barcodefile == "default:10x" ]]; then
        #restore logo in HTML template
        if [[ $verbose ]]; then
            echo "restore logo in summary HTML"
        fi
        if [[ -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html ]];then
            cp ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.html
        fi
        
        #restore cloupe generation
        if [[ $verbose ]]; then
            echo "restore cloupe" 
        fi
        ##list cloupe output as (not null)
        if [[ $verbose ]]; then
            echo "sed -i '/cloupe/s/null/CLOUPE_PREPROCESS\.output_for_cloupe/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i '/cloupe/s/null/CLOUPE_PREPROCESS\.output_for_cloupe/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro 
        ##add cloupe to outputs
        if [[ $verbose ]]; then
            echo "sed -i '/out cloupe *cloupe/ {s/^#*#//g}' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i '/out cloupe *cloupe/ {s/^#*#//g}' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        
        #restore 11 lines for cloupe preprocess call (all following steps are needed to be called together or the call will break)
        ##restore defining CLOUPE_PREPROCESS
        if [[ $verbose ]]; then
            echo "sed -i 's/^#*#@include "_cloupe_stages.mro"/@include "_cloupe_stages.mro"/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro"
        fi
        sed -i 's/^#*#@include "_cloupe_stages.mro"/@include "_cloupe_stages.mro"/g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro
        ##remove listing CLOUPE in output
        sed -i '/output_for_cloupe/s/^#*#//g' ${cellrangerpath}-cs/${cellrangerversion}/mro/*mro 
        ##remove calling CLOUPE_PREPROCESS
        ### iterate over all files calling CLOUPE_PREPROCESS
        for file in $(grep -l "call CLOUPE_PREPROCESS" ${cellrangerpath}-cs/${cellrangerversion}/mro/*.mro ); do
            #find start of CLOUPE_PREPROCESS call
            num=$(grep -n "call CLOUPE_PREPROCESS" $file | head -n 1 | cut -d":" -f 1)
            #find end of CLOUPE_PREPROCESS call
            num2=$(($num + $(tail -n $(echo $(($(wc -l $file | cut -f 1 -d " ") - $num))) $file | grep -n ")" | head -n 1 | cut -d":" -f 1)))
            if [[ $verbose ]]; then
               echo "lines ${num}-${num2} restored in $file"
            fi
            eval "sed -i '$(echo "${num},${num2}s/^#*#//g")' $file"
        done
    elif [[ $lastcall_p == "default:10x" ]] || [[ ! -f $lastcallfile ]]; then
        #remove logo from HTML template
        if [[ -f ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html ]];then
            #line of HTML in header is removed
            sed '/class="logo"/d' ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.backup.html > ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/webshim/template/navbar.html
        fi
        
        #disable cloupe generation
        if [[ $verbose ]]; then
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
        for file in $(grep -l "call CLOUPE_PREPROCESS" ${cellrangerpath}-cs/${cellrangerversion}/mro/*.mro ); do
            #find start of CLOUPE_PREPROCESS call
            num=$(grep -n "call CLOUPE_PREPROCESS" $file | head -n 1 | cut -d":" -f 1)
            #find end of CLOUPE_PREPROCESS call
            num2=$(($num +$(tail -n $(echo $(($(wc -l $file | cut -f 1 -d " ") - $num))) $file | grep -n ")" | head -n 1 | cut -d":" -f 1)))
            if [[ $verbose ]]; then
                echo "lines ${num}-${num2} commented out of $file"
            fi
            eval "sed -i '$(echo "${num},${num2}s/^/#/g")' $file"
        done
    fi
    
    #check for custom whitelist
    if [[ $custombarcodes == "true" ]]; then
        #disable detect chemistry check for custom whitelist
        if [[ $verbose ]]; then
            echo "disable detect chemistry check ..."
        fi
        sed -i "s/ raise NoChemistryFoundException/ return best_chem #raise NoChemistryFoundException/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/ (100\.0 \* best_frac/ #(100.0 * best_frac/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/return msg/return None #msg/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/return msg/return None #msg/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/check.py
    else
        if [[ $verbose ]]; then
            echo "restore detect chemistry check ..."
        fi
        
        #restore detect chemistry check for custom whitelist
        sed -i "s/ return best_chem \#raise NoChemistryFoundException/ raise NoChemistryFoundException/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/ \#(100\.0 \* best_frac/ (100.0 * best_frac/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/return None #msg/return msg/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/return None #msg/return msg/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/check.py
    fi
    
    #determine last barcode and UMI
    if [[ $lastcall_b == "" ]]; then
        old_bc_length=16
    else
        old_bc_length=$lastcall_b
    fi
    if [[ $lastcall_u == "" ]]; then
       old_umi_length=10
    else
       old_umi_length=$lastcall_u
    fi
    if [[ $verbose ]]; then
        echo "lastcall: b ${lastcall_b} u ${lastcall_u}; current: b ${barcodelength} u ${umilength}"
        echo "${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py"
    fi
    if [[ -z ${lastcall_b} ]]  || [[ -z ${lastcall_u} ]]; then
        old_rna_offset=26
    else
        old_rna_offset=`echo $((${lastcall_b} + ${lastcall_u}))` 2>%1 1> /dev/null
    fi
    new_rna_offset=`echo $((${barcodelength} + ${umilength}))`
    
    #convert barcodes back if last technology barcode greater than 16 bp
    if [[ $old_bc_length -gt 16 ]]; then
        if [[ $verbose ]]; then
           echo "barcode default restored to 16 bp"
        fi
        sed -i "s/'barcode_read_length': ${old_bc_length},/'barcode_read_length': 16,/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/'umi_read_offset': ${old_bc_length},/'umi_read_offset': 16,/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
    fi
    
    #convert UMI back if last technology UMI greater than 12 bp
    if [[ $old_umi_length -gt 12 ]]; then
        if [[ $verbose ]]; then
           echo "umi default restored to 10 bp"
        fi
        sed -i "s/'umi_read_length': ${old_umi_length},/'umi_read_length': 10,/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
    fi
    if [[ $old_rna_offset -gt 26 ]]; then
       if [[ $verbose ]]; then
           echo "RNA offset restored to 30 bp"
       fi
       sed -i "s/'rna_read_offset': ${old_rna_offset},/'rna_read_offset': 26,/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
       sed -i "s/'umi_read_length': ${old_umi_length},/'umi_read_length': 10,/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
    fi
    
    #convert barcodes if new technology greater than 16 bp
    if [[ $minlength -gt 16 ]]; then
        if [[ $verbose ]]; then
            echo "barcode length set to $minlength"
        fi
        sed -i "s/'barcode_read_length': 16,/'barcode_read_length': ${minlength},/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
        sed -i "s/'umi_read_offset': 16,/'umi_read_offset': ${minlength},/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
    fi
    
    #convert UMI back if new technology greater than 12 bp
    if [[ $umilength -gt 12 ]]; then
        if [[ $verbose ]]; then
            echo "umi length set to $umilength"
        fi
        sed -i "s/'umi_read_length': 10,/'umi_read_length': ${umilength},/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
    fi
    if [[ $new_rna_offset -gt 26 ]]; then
       if [[ $verbose ]]; then
           echo "RNA offset set to $new_rna_offset"
       fi
       sed -i "s/'rna_read_offset': 26,/'rna_read_offset': ${new_rna_offset},/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
       sed -i "s/'umi_read_length': 10,/'umi_read_length': ${umilength},/g" ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/chemistry.py
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
        if [[ $verbose ]]; then
            echo "  ... restoring 10x barcodes"
        fi
        #for version 2
        cp ${v2}.backup ${v2}
        #for version 3
        cp ${v3}.backup.gz ${v3}.gz
    else
        #for version 2
        cat ${barcodefile} > ${v2}
        echo "barcode adjust: $barcodeadjust"
        if [[ $barcodeadjust -gt 0 ]]; then
            sed -i "s/^.{${barcodeadjust}}//" ${v2} #Trim the first n characters from the beginning of the sequence and quality
        elif [[ 0 -gt $barcodeadjust ]]; then
            As=`printf '%0.sA' $(seq 1 $(($barcodeadjust * -1)))`
            echo As: $As
            sed -i "s/^/$As/" ${v2} #Trim the first n characters from the beginning of the quality
        fi
        #for version 3
        cat ${v2} > ${v3}
    fi
        if [[ -f ${v3} ]]; then
             gzip -f ${v3}
        fi
        if [[ -f translation/${v3} ]]; then
            gzip -f  translation/${v3}
        fi
        if [[ -f translation/${v3}.gz ]]; then
            rm translation/${v3}.gz
            zcat ${v3}.gz | awk -F , -v OFS="\t" '{print $1, "\t\t\t", $1}' > translation/${v3}
            gzip -f translation/${v3}
        fi
    echo " whitelist converted"
    
    echo "verbose $verbose"
    #change last call file
    if [[ $verbose ]]; then
        echo "setting last call as ..."
        echo "${barcodelength} ${umilength} ${barcodefile}"
    fi
    echo "${barcodelength} ${umilength} ${barcodefile}" > $lastcallfile
    
    cd - > /dev/null
    
    echo "setup complete"
fi
#########



#####checking job scheduler####
# check if cluster mode enabled
if [[ ${jobmode} != "local" ]]; then
    echo "running in cluster mode ${jobmode}"
    if [[ -d  $(dirname ${cellrangerpath})/martian-cs/*/jobmanagers ]]; then
        cd $(dirname ${cellrangerpath})/martian-cs/*/jobmanagers
        if [[ -f ${jobmode}.template ]]; then
            echo "cluster mode configured for ${jobmode}:" $(realpath ${jobmode}.template)
        else
            echo "attempting to use default template for cluster mode ${jobmode}"
            # copies default job template for SGE, Slurm, etc.
            if [[ -f ${jobmode}.template.example ]]; then
                cp -v ${jobmode}.template.example ${jobmode}.template
            fi
            echo "WARNING: if this does not work as expected, configure the scheduler by editing:"
            echo $(realpath ${jobmode}.template)
        fi
        #return to previus directory
        cd -
    fi
    if [[ -d  $(dirname ${cellrangerpath})/external/martian/jobmanagers ]]; then
        cd $(dirname ${cellrangerpath})/external/martian/jobmanagers
        if [[ -f ${jobmode}.template ]]; then
            echo "cluster mode configured for ${jobmode}:" $(realpath ${jobmode}.template)
        else
            echo "attempting to use default template for cluster mode ${jobmode}"
            # copies default job template for SGE, Slurm, etc.
            if [[ -f ${jobmode}.template.example ]]; then
                cp -v ${jobmode}.template.example ${jobmode}.template
            fi
            echo "WARNING: if this does not work as expected, configure the scheduler by editing:"
            echo $(realpath ${jobmode}.template)
        fi
        #return to previus directory
        cd -
    fi
else
    echo "running in local mode (no cluster configuration needed)"
fi
########



#####exiting when setup is all that is requested#####
if [[ $setup == "true" ]]; then
    lock=`cat $lockfile`
    if [[ $lock -ge 1 ]]; then
        lock=$(($lock - 1))
        echo $lock > $lockfile
    else
       rm -rf $lockfile
    fi
    echo " setup complete. exiting launch_universc.sh"
    exit 0
fi
##########



#####create directory with files fed to Cell Ranger#####
echo "creating a folder for all Cell Ranger input files ..."
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
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]] || [[ "$technology" == "vasa-drop" ]]; then
        #where converted "read1" is R2 in source files (corrected names for Cell Ranger)
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
    if [[ "$technology" == "indrop-v2" ]] || [[ "$technology" == "indrop-v3" ]] || [["$technology" == "vasa-drop" ]]; then
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
            #where converted "index1" is R2 in source files (corrected names for Cell Ranger)
            echo "using transcripts in Read 2 for ${technology}"
            to=`echo $to | sed -e "s/_R2_/_I1_/g"`
        fi
        
        if [[ $verbose ]]; then
            echo "$to"
        fi
        crI1s+=($to)
        
        echo "handling $fq ..."
        if [[ ! -f $to ]] || [[ $convert == true ]]; then
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
            #where converted "index2" is R3 in source files (corrected names for Cell Ranger)
            echo "using transcripts in Read 2 for ${technology}"
            to=`echo $to | sed -e "s/_R3_/_I2_/g"`
        fi
        
        if [[ $verbose ]]; then
            echo "$to"
        fi
        crI2s+=($to)
        
        echo "handling $fq ..."
        if [[ ! -f $to ]] || [[ $convert == true ]]; then
            if [[ $verbose ]]; then
                echo "cp -f $fq $to"
       	    fi
            cp -f $fq $to
        fi
    done
fi

if [[ ${#read3[@]} -ge 1 ]]; then
    crR3s=()
    if [[ $verbose ]]; then
        echo "Processing Index"
        echo "Fastqs: ${read3[@]}"
        echo "${read3[@]}"
    fi
    for fq in "${read3[@]}"; do
        if [[ $verbose ]]; then
            echo "$fq"
        fi
        to=`basename $fq`
        to="${crIN}/${to}"
        to=$(echo "$to" | sed 's/\.gz$//')
        
        if [[ $verbose ]]; then
            echo "$to"
        fi
        crR3s+=($to)
        
        echo "handling $fq ..."
        if [[ ! -f $to ]] || [[ $convert == true ]]; then
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
    echo "  barcodes: ${barcodeadjust} bp at its head"
    echo "  UMIs: ${umiadjust} bp at its tail" 
    
    echo " making technology-specific modifications ..."
    #CEL-Seq2 or VASA-plate: swap barcode and UMI
    ##https://github.com/BUStools/bustools/issues/4
    if [[ "$technology" == "celseq2" ]] || [[ "$technology" == "vasa-plate" ]]; then
        echo "  ... barcode and UMI swapped for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #swap UMI and barcode
            perl -pni -e "s/(.{6})(.{6})/\1\2/' if ($. % 2 == 0)" $convFile
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #BD Rhapsody: remove adapters
    if [[ "$technology" == "bd-rhapsody" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{9})ACTGGCCTGCGA(.{9})GGTAGCGGTGACA(.{9})(.{8})/ {
                s/.*(.{9})ACTGGCCTGCGA(.{9})GGTAGCGGTGACA(.{9})(.{8})/\1\2\3\4/g
                n
                n
                s/.*(.{9}).{12}(.{9}).{13}(.{9})(.{8})/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    #BD Rhapsody v2: remove adapters for enhanced beads (released in 2022)
    if [[ "$technology" == "bd-rhapsody" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*[ACG][CGT](.{9})GTGA(.{9})GACA(.{9})(.{8})/ {
                s/.*(.{9})GTGA(.{9})GACA(.{9})(.{8})/\1\2\3\4/g
                n
                n
                s/.*(.{9}).{4}(.{9}).{4}(.{9})(.{8})/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #C1, Quartz-Seq and RamDA-Seq: add mock UMI for non-UMI techniques
    if [[ "$technology" == "fluidigm-c1" ]] || [[ "$technology" == "bravo" ]] || [[ "$technology" == "c1-cage" ]] || [[ "$technology" == "ramda-seq" ]] || [[ "$technology" == "c1-ramda-seq" ]] || [[ "$technology" == "quartz-seq" ]]; then
        echo "  ... adding mock UMIs for ${technology}"
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            
            #detect barcode length from index sequence
            indexlength=$(($(head $index1([0]) -n 2 | tail -n 1 | wc -c) -1))
            barcodelength=$indexlength
            
            #detect whether index 2 files exist
            if [[ -f $(echo ${index2}([0])) ]]; then
                #detect barcode length from length of both index sequences
                index2length=$(($(head $index2([0]) -n 2 | tail -n 1 | wc -c) -1))
                barcodelength=$(($indexlength+$index2length))
                convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
                echo "  ... concatencate barcodes to R1 from I1 and I2 index files"
                #concatenate barcocdes from index to R1 as (bases 1-16 of the) barcode, moving (read to start at base 17-)
                perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir $crIN
            else
                barcodelength=$indexlength
                echo "  ... concatencate barcodes to R1 from I1 index files"
                #concatenate barcocdes from index to R1 as (bases 1-6 of the) barcode, moving (read to start at base 7-)
                perl ${CONCATENATEBARCODES} --additive ${convI1} --ref_fastq ${convR1} --out_dir $crIN
            fi
            
            #returns a combined R1 file with I1-R1 concatenated (I1 is cell barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
            
            if [[ $nonUMI == "true" ]]; then
                if [[ $verbose ]]; then
                    echo "adding mock UMI"
                fi
                #add mock UMI (count reads instead of UMI) barcodelength=6, umi_default=10
                perl ${ADDMOCKUMI} --fastq ${convR1} --out_dir $crIN --head_length $barcodelength --umi_length $umi_default
                mv $crIN/mock_UMI.fastq ${convR1}
                umilength=$umi_default
                umiadjust=0
                if [[ $chemistry == "SC3Pv3" ]]; then
                    chemistry="SC3Pv2"
                fi
                
                #avoid mock UMI being generated twice
                nonUMI=false
            fi
            
            if [[ "$technology" == "c1-cage" ]]; then
                #Add 10x TSO characters to the end of the sequence (removes 'NNNNNNNNTATAGGG')
                tsoS="TTTCTTATATGGG"
                tsoQ="IIIIIIIIIIIII"
                chemistry="SC5P-PE"
                if [[ $verbose ]]; then
                    echo "barcode: $barcodelength"
                    echo "umi: $umilength"
                fi
               
                #adjusting sequence data (removes 'NNNNNNNNTATAGGG')
                cmd=$(echo 'sed -E "2~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{8})TATAGGG/\1\2'$tsoS'/" '$convFile' > '${crIN}'/.temp')
                eval $cmd
                mv ${crIN}/.temp $convFile
                #Adjusting quality data (add n characters to the end of the quality)
                cmd=$(echo 'sed -E "4~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{8})(.{7})/\1\2'$tsoQ'/" '$convFile' > '${crIN}'/.temp')
                eval $cmd
                mv ${crIN}/.temp ${convFile}
            fi
            echo "  ${convFile} adjusted"
        done
    fi
    
    #ICELL8 version 2 (non-UMI technology)
    if [[ "$technology" == "icell8" ]] || [[ "$technology" == "icell8-5-prime" ]] || [[ "$technology" == "icell8-full-length" ]]; then
        echo "  ... filtering tagged reads for ${iechnology}"
        if [[ $verbose ]]; then
            echo "Note: ICELL8 v2 does not contain UMIs"
            echo "chemistry: $chemistry"
            echo "non-UMI: $nonUMI"
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            
            if [[ "$technology" == "icell8-full-length" ]]; then
                echo " ... remove internal reads for ${technology} by matching TSO sequence for UMI reads"
                #filter R1 reads by matching tag sequence AAGCAGTGGTATCAACGCAGAGTAC (same as SmartSeq2) and remove as an adapter 
                perl ${FILTERSMARTSEQREADUMI} --r1 ${convR1} --r2 ${convR2} --tag 'AAGCAGTGGTATCAACGCAGAGTAC' --out_dir $crIN
                echo "  ... trim tag sequence from R1"
                
                #returns R1 with tag sequence removed (left trim) starting with 8 bp UMI and corresponding reads for R2
                mv $crIN/parsed_R1.fastq ${convR1}
                mv $crIN/parsed_R2.fastq ${convR2}
                mv $crIN/parsed_I1.fastq ${convI1}
                mv $crIN/parsed_I2.fastq ${convI2}
            fi
            
            if [[ $nonUMI == "true" ]]; then
                if [[ $chemistry == "SC3P"* ]] && [[ "$technology" == "icell8" ]]; then
                    if [[ $verbose ]]; then
                        echo "remove UMI to replace with mock UMI"
                    fi
                    #remove inflated umi (to replace with mock UMI and count as reads)
                    sed -E '
                        /^@/ {
                        n
                        s/^(.{11})(.{14})(.*)/\1\3/g
                        n
                        n
                        s/^(.{11})(.{14})(.*)/\1\3/g
                        }' $convFile > ${crIN}/.temp
                    mv ${crIN}/.temp $convFile
                fi
                
                if [[ $verbose ]]; then
                    echo "adding mock UMI"
                fi
                #add mock UMI (count reads instead of UMI) barcodelength=10 or 11, umi_default=10 <- default for "icell8-5-prime" and "icell8-full-length"
                perl ${ADDMOCKUMI} --fastq ${convR1} --out_dir $crIN --head_length $barcodelength --umi_length $umi_default
                mv $crIN/mock_UMI.fastq ${convR1}
                umilength=$umi_default
                umiadjust=0
                if [[ $chemistry == "SC3Pv3" ]]; then
                    chemistry="SC3Pv2"
                fi
                
                #avoid mock UMI being generated twice
                nonUMI=false
            fi
            
            #subroutine for icell8-5-prime (5' scRNA with TSO adapters) or icell8-full-length (assumes barcodes and mock UMI added above)
            if [[ $chemistry == "SC5P"* ]] || [[ $chemistry == "five"* ]]; then
                #remove TSO adapters (from ends)
                sed -E '
                    /^(.{10})(.{10})GCGTCGTGTAGGG/ {
                    s/^(.{10})(.{10})GCGTCGTGTAGGG/\1\2GGG/g
                    n
                    n
                    s/^(.{10})(.{10})(.{10})(.{3})/\1\2\4/g
                }'  $convFile > ${crIN}/.temp
                mv ${crIN}/.temp $convFile
                sed -E '
                    /^(.{10})(.{10})TGTGAGAAAGGG/ {
                    s/^(.{10})(.{10})TGTGAGAAAGGG/\1\2GGG/g
                    n
                    n
                    s/(.{10})(.{10})(.{9})(.{3})/\1\2\4/g
                }'  $convFile > ${crIN}/.temp
                mv ${crIN}/.temp $convFile
                #convert TSO to expected length for 10x 5' (TSS in R1 from base 39)
                echo " handling $convFile ..."
                tsoS="TTTCTTATATGGG" #<- confirm tag sequence with Takara reps
                tsoQ="IIIIIIIIIIIII"
                #Add 10x TSO characters to the end of the sequence (remove 'GGG' linker)
                cmd=$(echo 'sed -E "2~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{3})/\1\2'$tsoS'/" '$convFile' > '${crIN}'/.temp') #<- confirm tag sequence (GGG) with Takara reps
                if [[ $verbose ]]; then
                    echo technology $technology
                    echo barcode: $barcodelength
                    echo umi: $umilength
                    echo $cmd
                fi
                #run command with barcode and umi length, e.g.,: sed -E "2~4s/(.{10})(.{10})(.{3})(.*)/\1\2$tsoS\4/" $convFile > ${crIN}/.temp
                eval $cmd
                mv ${crIN}/.temp $convFile
                #Add n characters to the end of the quality
                cmd=$(echo 'sed -E "4~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{3})/\1\2'$tsoQ'/" '$convFile' > '${crIN}'/.temp')
                #run command with barcode and umi length, e.g.,: sed -E "4~4s/(.{10})(.{10})(.{3})(.*)/\1\2$tsoQ\4/" $convFile > ${crIN}/.temp
                eval $cmd
                mv ${crIN}/.temp $convFile
            fi
            echo "  ${convFile} adjusted"
       done
    fi
    
    #inDrops: remove adapter (see links below for details)
    ##https://github.com/BUStools/bustools/issues/4
    ##https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html
    ##https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup
    #note that adapters do not have to be removed for the dual-indexed inDrops-v3
    if [[ "$technology" == "indrop-v1" ]] || [[ "$technology" == "indrop-v2" ]]; then
        echo "  ... remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapter if present
            sed -E '
                /^.*(.{8})GA.................CTT(.{8})(.{6}).*$/ {
                s/^.*(.{8})GA.................CTT(.{8})(.{6}).*$/\1\2\3/g
                n
                n
                s/^(.{8}).{22}(.{8})(.{6}).*/\1\2\3/g
                }' $convFile > ${crIN}/.temp
            #remove linker between barcode and UMI
            echo "  ... barcode and UMI linker removed for ${technology}"
            mv ${crIN}/.temp $convFile
        done
    fi
    #inDrops: migrate dual indexes to barcode
    #https://github.com/BUStools/bustools/issues/4
    if [[ "$technology" == "indrop-v3" ]]; then
        echo "  ... processsing for ${technology}"
         if [[ $verbose ]]; then
             echo "Note: inDrops v3 should be demultiplex by sample index I2 (R3 or 4) if multiple samples are sequenced"
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R2/$1_R1/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R2/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R2/$1_I2/' )
            
            #(R1 -> R2; R2 -> I1; R3 -> I2; R4 -> R1)
            #v3: summer 2016 redesign requiring manual demultiplexing.
            #R1 is the biological read (R2).
            #R2 (i7) carries the first half of the gel barcode (I1). <- which cell
            #R3 (i5) carries the library index (I2). <- which sample/organism
            #R4 the second half of the gel barcode, the UMI and a fraction of the polyA tail (R1).
            #returns R1 with tag sequence removed (left trim) starting with 8 bp UMI and corresponding reads for I1, I2, and R2
            
            echo "  ... concatencate barcodes to R1 from I1 index file"
            #concatenate barcocdes from dual indexes to R1 as (bases 1-8 of the) barcode (bases 1-16), moving UMI to (17-22)
            #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as an adapters
            perl ${CONCATENATEBARCODES} --additive ${convI1} --ref_fastq ${convR1} --out_dir $crIN
            
            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
        done
    fi
    
    #Microwell-Seq: remove linkers
    if [[ "$technology" == "microwellseq" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{6})CGACTCACTACAGGG(.{6})TCGGTGACACGATCG(.{6})(.{6})/ {
                s/.*(.{6})CGACTCACTACAGGG(.{6})TCGGTGACACGATCG(.{6})(.{6})/\1\2\3\4/g
                n
                n
                s/.*(.{6}).{15}(.{6}).{15}(.{6})(.{6})$/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #PIP-Seq: remove adapter and correct phase blocks
    if [[ "$technology" == "pip-seq-v0" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{8})AACC(.{8})ACAG(.{8})CCTATTCGAG(.{8})/ {
                s/.*(.{8})AACC(.{8})ACAG(.{8})ACG(.{8})CCTATTCGAG(.{8})/\1\2\3\4/g
                n
                n
                s/.*(.{8}).{4}(.{8}).{4}(.{8}).{10}(.{8})/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    if [[ "$technology" == "pip-seq-v1" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{8})GAGTGATTGCTTGTGACGCCTT(.{8})(.{6})/ {
                s/.*(.{8})GAGTGATTGCTTGTGACGCCTT(.{8})(.{6})/\1\2\3/g
                n
                n
                s/.*(.{8}).{22}(.{8})(.{6})/\1\2\3/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    if [[ "$technology" == "pip-seq-v2" ]]; then
        echo "  ... remove adapters for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /^(.{8})ATGCATC(.{8})CCTCGAG(.{8})(.{12})/ {
                s/^(.{8})ATGCATC(.{8})CCTCGAG(.{8})(.{12})/\1\2\3\4/g
                n
                n
                s/^(.{8}).{7}(.{8}).{7}(.{8})(.{12})/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    if [[ "$technology" == "pip-seq-v3" ]]; then
        echo "  ... remove adapters for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /^(.{8})ATG(.{8})GAG(.{8})TCGAG(.{8})(.{12})/ {
                s/^(.{8})ATG(.{8})GAG(.{8})TCGAG(.{8})(.{12})/\1\2\3\4\5/g
                n
                n
                s/^(.{8}).{3}(.{8}).{3}(.{8}).{5}(.{8})(.{12})/\1\2\3\4\5/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    if [[ "$technology" == "pip-seq-v4" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{8})ATG(.{8})GAG(.{8})TCGAG(.{8})(.{12})/ {
                s/.*(.{8})ATG(.{8})GAG(.{8})TCGAG(.{8})(.{12})/\1\2\3\4\5/g
                n
                n
                s/.*(.{8}).{3}(.{8}).{3}(.{8}).{5}(.{8})(.{12})/\1\2\3\4\5/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi


    #Quartz-Seq2: remove adapter
    if [[ "$technology" == "quartz-seq2-384" ]]; then
        for convFile in "${convFiles[@]}"; do
        echo "  ... remove adapter for ${technology}"
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
        echo "  ... remove adapter for ${technology}"
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
    if [[ "$technology" == "sciseq2" ]]; then
        echo "  ... remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapter if detected (two-level indexing)
             ##TruSeq adapter: ACGACGCTCTTCCGATCT
            sed -E '
                /^ACGACGCTCTTCCGATCT/ {
                s/^ACGACGCTCTTCCGATCT(.{18})/\1/g
                n
                n
                s/^(.{18})(.{18})/\2/g
                }'  $convFile > ${crIN}/.temp
             mv ${crIN}/.temp $convFile
             #swap barcode and UMI
            echo "  ... barcode and UMI swapped for ${technology}"
            perl -pni -e "s/(.{8})(.{10})/\1\2/' if ($. % 2 == 0)" $convFile
            
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            
            echo "  ... concatencate barcodes to R1 from I1 and I2 index files"
            #concatenate barcocdes from dual indexes to R1 as (bases 1-20 of the) barcode, moving RT barcode (21-30) UMI to (31-38)
            #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as an adapters
            perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}

            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
        done
    fi
    
    
    if [[ "$technology" == "sciseq3" ]]; then
        echo "  ... remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapter if detected (and keep hairpin/tn5 barcode)
            ##TruSeq adapter: ACGACGCTCTTCCGATCT
            sed -E '
                /^ACGACGCTCTTCCGATCT/ {
                s/^ACGACGCTCTTCCGATCT//g
                n
                n
                s/^.{18}//g
                }' $convFile > ${crIN}/.temp
           mv ${crIN}/.temp $convFile
           #remove linker (10 bp barcodes)
           cat $convFile | sed -E '
                /^(.{9})CAGAGC/ {
                s/^(.{9})CAGAGC(.{18})/T\1CAGAGC\2/g
                n
                n
                s/^(.{9})(.{6})(.{18})/F\1\2\3/g
                }' |
           #remove linker (9 bp barcodes)
           sed -E '
                /^(.{10})CAGAGC/ {
                s/^(.{10})CAGAGC/\1/g
                n
                n
                s/^(.{10})(.{6})/\1/g
                }'  > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
            #swap barcode and UMI
            echo "  ... barcode and UMI swapped for ${technology}"
            perl -pni -e "s/(.{10})(.{8})(.{10})/\1\3\2/' if ($. % 2 == 0)" $convFile
            
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            
            echo "  ... concatencate barcodes to R1 from I1 and I2 index files"
            #concatenate barcocdes from dual indexes to R1 as (bases 1-20 of the) barcode, moving RT barcode (21-30) UMI to (31-38)
            #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as an adapters
            perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}
            
            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
        done
    fi
    
    if [[ "$technology" == "scifiseq" ]]; then
        echo "  ... remove adapter for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove adapter if detected (and keep hairpin/tn5 barcode)
            ##TruSeq adapter: ACGACGCTCTTCCGATCT
            sed -E '
                /^ACACTCTTTCCCTACACGACGCTCTTCCGATCT/ {
                s/^ACACTCTTTCCCTACACGACGCTCTTCCGATCT//g
                n
                n
                s/^.{33}//g
                }' $convFile > ${crIN}/.temp
           mv ${crIN}/.temp $convFile
           #remove linker and swap RT barcode and UMI (removes one base on either end of barcode): 11 bp barcode not 13 bp
           echo "  ... barcode and UMI swapped for ${technology}"
           sed -E '
                /^(.{8}).(.{11)./ {
                s/^(.{8}).(.{11)./\2\1/g
                n
                n
                s/^(.{8}).(.{11)./\2\1/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
            
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            
            echo "  ... concatencate 10x ATAC barcodes to R1 from I2 index files"
            #concatenate barcocdes from dual indexes to R1 as (bases 1-16 of the 27 bp) barcode, moving RT barcode (17-27) UMI to (28-35)
            #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as an adapters
            perl ${CONCATENATEBARCODES} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}
            
            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
        done
    fi
    
    #STRT-Seq
    if [[ "$technology" == "strt-seq" ]] || [[ "$technology" == "strt-seq-c1" ]] || [[ "$technology" == "strt-seq-2i" ]]; then
        echo "  ... processsing for ${technology}"
        if [[ $verbose ]]; then
            if [[ "$technology" == "strt-seq" ]]; then
                echo "Note: STRT-Seq does not contain UMIs"
            fi
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            
            if [[ "$technology" == "strt-seq" ]]; then
                # add mock UMI (count reads instead of UMI) barcodelength=6, umi_default=10
                echo "  ... generate mock UMI for compatibility"
                perl ${ADDMOCKUMI} --fastq ${convR1} --out_dir ${crIN} --head_length ${barcodelength} --umi_length ${umi_default}
                umilength=${umi_default}
                umiadjust=0
                chemistry="SC5P-R1"
                
                #returns a combined R1 file with barcode and mock UMI
                ##6 bp barcode, 10 bp UMI, GGG for TSO
                mv $crIN/mock_UMI.fastq ${convR1}
            fi
            if [[ "$technology" == "strt-seq-c1" ]]; then
                convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
                chemistry="SC5P-R1"
                echo "  ... concatencate barcodes to R1 from I1 index files"
                #concatenate barcocdes from I1 index to R1 as barcode (bases 1-8)
                perl ${CONCATENATEBARCODES} --additive ${convI1} --ref_fastq ${convR1} --out_dir ${crIN}
            fi
            if [[ "$technology" == "strt-seq-2i" ]]; then
                convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
                convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
                chemistry="SC5P-R1"
                echo "  ... concatencate barcodes to R1 from I1 and I2index files"
                #concatenate barcocdes from I1 index to R1 as barcode (bases 1-8)
                perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}
            fi
            
            #convert TSO to expected length for 10x 5' (TSS in R1 from base 39)
            echo " handling ${convFile} ..."
            tsoS="TTTCTTATATGGG"
            tsoQ="IIIIIIIIIIIII"
            #Add 10x TSO characters to the end of the sequence
            cmd=$(echo 'sed -E "2~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{3})/\1\2'$tsoS'/" '$convFile' > '${crIN}'/.temp')
            if [[ $verbose ]]; then
                echo technology $technology
                echo barcode: $barcodelength
                echo umi: $umilength
                echo $cmd
            fi
            #run command with barcode and umi length, e.g.,: sed -E "2~4s/(.{16})(.{8})(.{3})(.*)/\1\2$tsoS\4/" $convFile > ${crIN}/.temp
            eval $cmd
            mv ${crIN}/.temp $convFile
            #Add n characters to the end of the quality
            cmd=$(echo 'sed -E "4~4s/(.{'$barcodelength'})(.{'${umilength}'})(.{3})/\1\2'$tsoQ'/" '$convFile' > '${crIN}'/.temp')
            # run command with barcode and umi length, e.g.,: sed -E "4~4s/(.{16})(.{8})(.{3})(.*)/\1\2$tsoQ\4/" $convFile > ${crIN}/.temp
            eval $cmd
            mv ${crIN}/.temp $convFile
            echo "  ${convFile} adjusted"
        done
    fi

    #SureCell: remove adapter and correct phase blocks
    ##https://github.com/Hoohm/dropSeqPipe/issues/42
    if [[ "$technology" == "surecell" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC/ {
                s/.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC/\1\2\3\4/g
                n
                n
                s/.*(.{6}).{15}(.{6}).{15}(.{6}).{3}(.{8}).{3}/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    
    #SPLiT-Seq: correct phase blocks and swap barcode and UMI (if a whitelist and 24 bp barcode can be supported)
    ##https://github.com/hms-dbmi/dropEst/issues/80
    ##https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup
    if [[ "$technology" == "splitseq" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /^(.{8})CGAATGC........CTCAAGCACGTG[AG]AT(.{8})AGTC[GT]T[AC]..........[AC]C.CA[GT]C..[AC]C..(.{8})(.{8}).*/ {
                s/^(.{8})CGAATGC........CTCAAGCACGTG[AG]AT(.{8})AGTC[GT]T[AC]..........[AC]C.CA[GT]C..[AC]C..(.{8})(.{8}).*/\1\2\3\4/g
                n
                n
                s/^(.{8}).{30}(.{8}).{30}(.{8})(.{8}).*/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
            #remove phase blocks and linkers (reverse complement if R2 matched)
            sed -E '
                /[ATCGN][ATGCGN](.{8})(.{8}).{30}(.{8}).{30}(.{8}).*/ {
                s/..(.{8})(.{8}).{30}(.{8}).{30.(.{8}).*/\4\3\2\1/g
                n
                n
                s/..(.{8})(.{8}).{30}(.{8}).{30}(.{8}).*/\4\3\2\1/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi
    ## accounts for shorter adapters for SLiPT-Seq v2.1
    if [[ "$technology" == "splitseq2" ]]; then
        echo "  ... remove adapter and phase blocks for ${technology}"
        for convFile in "${convFiles[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /^(.{8})C[GC]CC...CTCCCGCCCGTG[CG]CT(.{8})AGTC[GT]T[AC]..........[AC]C.CA[GT]C..[AC]C..(.{8})(.{10}).*/ {
                s/^(.{8})C[GC]CC...CTCCCGCCCGTG[CG]CT(.{8})AGTC[GT]T[AC]..........[AC]C.CA[GT]C..[AC]C..(.{8})(.{10}).*/\1\2\3\4/g
                n
                n
                s/^(.{8}).{24}(.{8}).{30}(.{8})(.{10}.*$/\1\2\3\4/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
            #remove phase blocks and linkers (reverse complement if R2 matched)
            sed -E '
                /^(.{10})(.{8})..G[GT]..G[AC]TG.G[GT]..........[GT]A[AC]GACT(.{8})AT[CT]CACGTGCTTGAG...GT[GC]G(.{8}).*/ {
                s/^(.{10})(.{8})..G[GT]..G[AC]TG.G[GT]..........[GT]A[AC]GACT(.{8})AT[CT]CACGTGCTTGAG...GT[GC]G(.{8}).*/\4\3\2\1/g
                n
                n
                s/^(.{10})(.{8}).{30}(.{8}).{24}(.{8}).*/\4\3\2\1/g
                }' $convFile > ${crIN}/.temp
            mv ${crIN}/.temp $convFile
        done
    fi

    #SCRB-Seq: remove adapter
    ##https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html
    if [[ "$technology" == "scrbseq" ]]; then
        echo "  ... remove adapter for ${technology}"
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
    
    #Smart-Seq2 (or Smart-Seq)
    if [[ "$technology" == "smartseq2" ]];then
        echo "  ... processsing for ${technology}"
        if [[ $verbose ]]; then
            echo "Note: SmartSeq-2 does not contain UMIs"
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            
            #remove R1 adapter if detected
            sed -E '
                /^AGATGTGTATAAGAGACAG/ {
                s/^(.{19})//g
                n
                n
                s/^(.{19})//g
                }' $convR1 > ${crIN}/.temp
            mv ${crIN}/.temp $convR1
            #remove R2 adapter if detected
            sed -E '
                /CTGTCTCTTATACACATCT$/ {
                s/(.{19})$//g
                n
                n
                s/(.{19})$//g
                }' $convR2 > ${crIN}/.temp
            mv ${crIN}/.temp $convR2
            
            if [[ ${chemistry} == "SC5P-R1" ]]; then
                echo " ... remove internal reads for ${technology} by matching TSO sequence for UMI reads"
                #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as an adapters 
                perl ${FILTERSMARTSEQREADUMI} --r1 ${convR1} --r2 ${convR2} --i1 ${convI1} --i2 ${convI2} --tag 'AAGCAGTGGTATCAACGCAGAGTACGG' --out_dir ${crIN}
                echo "  ... trim tag sequence from R1"

                # returns R1 with tag sequence removed (left trim) starting with a 13 bp TSO and corresponding reads for I1, I2, and R2
                mv $crIN/parsed_R1.fastq ${convR1}
                mv $crIN/parsed_R2.fastq ${convR2}
                mv $crIN/parsed_I1.fastq ${convI1}
                mv $crIN/parsed_I2.fastq ${convI2}
                if [[ $verbose ]]; then
                    cp ${convR1} $crIN/parsed_R1.fastq
                    cp ${convR2} $crIN/parsed_R2.fastq
                    cp ${convI1} $crIN/parsed_I1.fastq
                    cp ${convI2} $crIN/parsed_I2.fastq
                fi
            else
                #if [[ ${chemistry} == "SC3Pv2" ]] || [[ ${chemistry} == "SC5P-PE" ]];
                # remove tag sequence adapter (first occurence only)
                sed -E '
                     /^AAGCAGTGGTATCAACGCAGAGTACGG/ {
                     s/^AAGCAGTGGTATCAACGCAGAGTACGG//g
                     n
                     n
                     s/^.{27}//g
                     }' $convFile > ${crIN}/.temp
                # returns R1 with tag sequence removed  
                mv ${crIN}/.temp $convFile
            fi
            
            echo "  ... concatencate barcodes to R1 from I1 and I2 index files"
            #concatenate barcocdes from dual indexes to R1 as barcode (bases 1-16)
            perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}
            
            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            ##16 bp barcode, GGG for TSO, no UMI
            mv $crIN/Concatenated_File.fastq ${convR1}
            
            if [[ $nonUMI == "true" ]]; then
                echo "  ... add mock UMI to count reads for to R1 files"
                #add mock UMI (count reads instead of UMI) barcodelength=16, umi_default=10
                perl ${ADDMOCKUMI} --fastq ${convR1} --out_dir ${crIN} --head_length ${barcodelength} --umi_length ${umi_default}
                umilength=${umi_default}
                umiadjust=0
            
                # returns a combined R1 file with barcode and mock UMI
                ##16 bp barcode, 10 bp UMI, GGG for TSO
                mv $crIN/mock_UMI.fastq ${convR1}
            fi

            # skip adding TSO for 3' chemistry (and 5' R1 which includes this in tag filtering)
            if [[ ${chemistry} != "SC3Pv2" ]] && [[ ${chemistry} != "SC3Pv3" ]] && [[ ${chemistry} != "SC5P-R1" ]] && [[ ${chemistry} != "auto" ]]; then
                #convert TSO to expected length for 10x 5' (TSS in R1 from base 39)
                echo " handling $convFile ..."
                tsoS="TTTCTTATATGGG"
                tsoQ="IIIIIIIIIIIII"
                #Add 10x TSO characters to the end of the sequence
                cmd=$(echo 'sed -E "2~4s/(.{'$barcodelength'})(.{'${umilength}'})/\1\2'$tsoS'/" '$convFile' > '${crIN}'/.temp')

                if [[ $verbose ]]; then
                    echo technology $technology
                    echo barcode: $barcodelength
                    echo umi: $umilength
                    echo $cmd
                fi

                #run command with barcode and umi length, e.g.,: sed -E "2~4s/(.{16})(.{8})(.*)/\1\2$tsoS\4/" $convFile > ${crIN}/.temp
                eval $cmd
                mv ${crIN}/.temp $convFile
                #Add n characters to the end of the quality
                cmd=$(echo 'sed -E "4~4s/(.{'$barcodelength'})(.{'${umilength}'})/\1\2'$tsoQ'/" '$convFile' > '${crIN}'/.temp')
                #run command with barcode and umi length, e.g.,: sed -E "4~4s/(.{16})(.{8})(.*)/\1\2$tsoQ\4/" $convFile > ${crIN}/.temp
                eval $cmd
                mv ${crIN}/.temp $convFile
            fi
            echo "  ${convFile} adjusted"
        done
    fi
    
    #Smart-Seq3
    if [[ "$technology" == "smartseq3" ]];then
        echo "  ... processsing for ${technology}"
        if [[ $verbose ]]; then
            echo "Note: SmartSeq-3 requires additional filtering for UMIs"
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            
            #filter UMI reads by matching tag sequence ATTGCGCAATG (bases 1-11 of R1) and remove as adapters 
            echo "  ... parsing R1 reads with tag sequence and inserting 10x TSO"
            perl ${FILTERSMARTSEQREADUMI} --r1 ${convR1} --r2 ${convR2} --i1 ${convI1} --i2 ${convI2} --tag 'ATTGCGCAATG' --out_dir ${crIN}
             
            #returns R1 with tag sequence removed (left trim) starting with 8bp UMI and 13 bp TSO and corresponding reads for I1, I2, and R2
            mv $crIN/parsed_R1.fastq ${convR1}
            mv $crIN/parsed_R2.fastq ${convR2}
            mv $crIN/parsed_I1.fastq ${convI1}
            mv $crIN/parsed_I2.fastq ${convI2}
            if [[ $verbose ]]; then
                cp ${convR1} $crIN/parsed_R1.fastq
                cp ${convR2} $crIN/parsed_R2.fastq
                cp ${convI1} $crIN/parsed_I1.fastq
                cp ${convI2} $crIN/parsed_I2.fastq
            fi
            
            #concatenate barcocdes from dual indexes to R1 as barcode (bases 1-16)
            echo "  ... concatencate barcodes to R1 from I1 and I2 index files"
            perl ${CONCATENATEBARCODES} --additive ${convI1} --additive ${convI2} --ref_fastq ${convR1} --out_dir ${crIN}
            
            #returns a combined R1 file with I1-I2-R1 concatenated (I1 and I2 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
            if [[ $verbose ]]; then
                cp ${convR1} $crIN/Concatenated_File.fastq
            fi
            echo "  ${convFile} adjusted"
        done
    fi

    #VASA-drop
    if [[ "$technology" == "vasa-drop" ]];then
        echo "  ... processsing for ${technology}"
        if [[ $verbose ]]; then
            echo "Note: VASA-drop contains reads barcodes in I1 and R2 (switched to R1 already)"
        fi
        for convFile in "${convFiles[@]}"; do
            read=$convFile
            convR1=$read
            convR2=$(echo $read | perl -pne 's/(.*)_R1/$1_R2/' )
            convI1=$(echo $read | perl -pne 's/(.*)_R1/$1_I1/' )
            convI2=$(echo $read | perl -pne 's/(.*)_R1/$1_I2/' )
            convI1bcs=$(echo $read | perl -pne 's/(.*)_R1/$1_I1_bcs/' )

            #create a temporary files with I1 barcodes extracted
            cut -c 1-8 ${convI1} > $crIN/parsed_I1.fastq ${convI1bcs}

            #concatenate barcocdes from index I1 to R1 as barcode (bases 1-8 from I2; bases 9-16 from R2)
            echo "  ... concatencate barcodes to R1 from I1 index files"
            perl ${CONCATENATEBARCODES} --additive ${convI1bcs} --ref_fastq ${convR1} --out_dir ${crIN}

            #returns a combined R1 file with I1-R1 concatenated (I1 are R1 barcode)
            mv $crIN/Concatenated_File.fastq ${convR1}
            if [[ $verbose ]]; then
                cp ${convR1} $crIN/Concatenated_File.fastq
                cp ${convI1bcs} $crIN/I1_barcodes.fastq
            fi
            rm -rf ${convI1bcs}
            echo "  ${convFile} adjusted"
        done
    fi
    
    #converting barcodes
    echo " adjusting barcodes of R1 files"
    if [[ $barcodeadjust != 0 ]]; then
        if [[ $barcodeadjust -gt 0 ]]; then
            for convFile in "${convFiles[@]}"; do
                echo " handling $convFile ..."
                perl -pni -e "s/^.{${barcodeadjust}}//g if ($. % 2 == 0)" $convFile #Trim the first n characters from the beginning of the sequence and quality
                echo "  ${convFile} adjusted"
            done
        elif [[ 0 -gt $barcodeadjust ]]; then
            for convFile in "${convFiles[@]}"; do
                echo " handling $convFile ..."
                toS=`printf '%0.sA' $(seq 1 $(($barcodeadjust * - 1)))`
                toQ=`printf '%0.sI' $(seq 1 $(($barcodeadjust * - 1)))`
                perl -pni -e "s/^/$toS/g if ($. % 4 == 2)" $convFile #Trim the first n characters from the beginning of the sequence
                perl -pni -e "s/^/$toQ/g if ($. % 4 == 0)" $convFile #Trim the first n characters from the beginning of the quality
                echo "  ${convFile} adjusted"
            done
        fi
    fi
    
    #replace UMI with mock UMI to count reads (for technologies not already containing mock UMI)
    if [[ $technology != "icell8" ]] && [[ $technology != "ramda-seq" ]] && [[ $technology != "quartz-seq" ]] && [[ $technology != "smartseq" ]] && [[ $technology != "smartseq2" ]] && [[ $technology != "strt-seq" ]]; then
        if [[ $nonUMI == "true" ]]; then
            echo "***WARNING: Removing true UMI and replacing with Mock UMI. This is not recommended unless integrating with non-UMI data***"
            echo "NOTICE: results will show read counts and not UMI counts"
             
            for convFile in "${convFiles[@]}"; do
                convR1=$convFile
                #remove inflated umi (to replace with mock and count as reads) by non-standard evaluation to depend on variable umi-length
                cmd=$(echo 'sed -E "2~4s/\^(.{'$barcodelength'})(.{'${umilength}'})(.\*)$/\1\3/" '$convFile' > '${crIN}'/.temp')
                if [[ $verbose ]]; then
                    echo technology $technology
                    echo barcode: $barcodelength
                    echo umi: $umilength
                    echo $cmd
                fi
                eval $cmd
                mv ${crIN}/.temp $convFile
                
                if [[ $chemistry == "SC3Pv3" ]]; then
                     chemistry="SC3Pv2"
                fi
                perl ${ADDMOCKUMI} --fastq ${convR1} --out_dir $crIN --head_length ${barcodelength} --umi_length ${umi_default}
                umilength=$umi_default
                umiadjust=0
                if [[ $chemistry == "SC3Pv3" ]]; then
                    chemistry="SC3Pv2"
                fi
                
                #returns a combined R1 file with barcode and mock UMI
                ## barcode, 10 bp UMI, followed by TSO (if applicable)
                mv $crIN/mock_UMI.fastq ${convR1}
            done
        fi
    fi
    
    #convert UMI
    echo " adjusting UMIs of R1 files"
    # check if original UMI is shorter than default
    if [[ 0 -gt $umiadjust ]]; then
        for convFile in "${convFiles[@]}"; do
            echo " handling $convFile ..."
            toS=`printf '%0.sA' $(seq 1 $(($umiadjust * - 1)))`
            toQ=`printf '%0.sI' $(seq 1 $(($umiadjust * - 1)))`
            #compute length of adjusted barcode + original UMI
            keeplength=`echo $((${barcode_default} + ${umi_default} - ($umiadjust * - 1)))`
            #Add n characters to the end of the sequence
            perl -pni -e "s/(.{$keeplength})(.*)/\1$toS\2/ if ($. % 4 == 2)" $convFile
            #Add n characters to the end of the quality
            perl -pni -e "s/(.{$keeplength})(.*)/\1$toQ\2/ if ($. % 4 == 0)" $convFile
            echo "  ${convFile} adjusted"
        done
    fi
    
    # check if original UMI is longer than default
    if [[ 0 -lt $umiadjust ]]; then
        for convFile in "${convFiles[@]}"; do
            echo " handling $convFile ..."
            #compute length of adjusted barcode + original UMI
            targetlength=`echo $((${barcode_default} + ${umi_default}))`
            #Remove n characters to the end of the sequence and quality score
            perl -pni -e "s/(.{$targetlength})(.{$umiadjust})(.*)/\1\3/ if ($. % 2 == 0)" $convFile
            mv ${crIN}/.temp $convFile
            echo "  ${convFile} adjusted"
        done
    fi
fi
##########



#####setting parameters for Cell Ranger#####
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

if [[ $verbose ]]; then
    echo $chemistry
fi

if [[ $chemistry == "SC5P"* ]] || [[ $chemistry == "five"* ]] || [[ $chemistry == "auto" ]]; then
    r=""
else
    r="--r1-length=$totallength"
fi

#outputting command
echo "running Cell Ranger ..."
echo ""
echo "#####Cell Ranger command#####"

start=`date +%s`
echo "cellranger count --id=$id\\
        --fastqs=$crIN\\
        --lanes=$LANE\\
        $r\\
        --chemistry=$chemistry\\
        --transcriptome=$reference\\
        --sample=$SAMPLE\\
        $d\\
        $n\\
        $j\\
        $l\\
        $m\\
        $bam\\
        $secondary\\
        $library\\
        $introns\\
        $noexit\\
        $nopreflight
"
echo "##########"
##########



#####running Cell Ranger#####
cellranger count --id=$id \
        --fastqs=$crIN \
        --lanes=$LANE \
        $r \
        --chemistry=$chemistry \
        --transcriptome=$reference \
        --sample=$SAMPLE \
        $d \
        $n \
        $j \
        $l \
        $m \
        $bam \
        $secondary \
        $library \
        $introns \
        $noexit \
        $nopreflight
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



#####remove files if universc is not running elsewhere#####
echo "updating .lock file"

#remove current job from counter (successfully completed)
lock=`cat ${cellrangerpath}-cs/${cellrangerversion}/lib/python/cellranger/barcodes/.lock`
lock=$(($lock - 1))
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

