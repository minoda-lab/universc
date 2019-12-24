#!/bin/bash

install=false

#####cellrenger version#####
cellrangerpass=`which cellranger`
if [[ -z $cellrangerpass ]]; then
    echo "cellranger command is not found."
    exit 1
fi
ver_info=`paste -d "\n" <(cellranger count --version) <(echo conversion script version 0.2.2) | head -n 3 | tail -n 2`
##########



#####locate launch_universc.sh for importing barcodes######
SOURCE="${BASH_SOURCE[0]}"
while [[ -h "$SOURCE" ]]; do #resolve $SOURCE until the file is no longer a symlink
    TARGET="$(readlink "$SOURCE")"
    if [[ $TARGET == /* ]]; then
        echo "SOURCE '$SOURCE' is an absolute symlink to '$TARGET'"
        SOURCE="$TARGET"
    else
        SCRIPT_DIR="$( dirname "$SOURCE" )"
        echo "SOURCE '$SOURCE' is a relative symlink to '$TARGET' (relative to '$SCRIPT_DIR')"
        SOURCE="$SCRIPT_DIR/$TARGET" #if $SOURCE is a relative symlink, we need to resolve it relative to the path where the symlink file was located
    fi
done
RDIR="$( dirname "$SOURCE" )"
SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
if [[ $RDIR != $SCRIPT_DIR ]]; then
    echo "DIR '$RDIR' resolves to '$SCRIPT_DIR'"
fi
echo "Running launch_universc.sh in '$SCRIPT_DIR'"
##########



#####usage statement#####
help='
Usage:
  bash '$(basename $0)' -R1 FILE1 -R2 FILE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -R1 READ1_LANE1 READ1_LANE2 -R2 READ2_LANE1 READ2_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -f SAMPLE_LANE -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -f SAMPLE_LANE1 SAMPLE_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash '$(basename $0)' -v
  bash '$(basename $0)' -h
  bash '$(basename $0)' -t TECHNOLOGY --setup

Convert sequencing data (FASTQ) from Nadia or iCELL8 platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
  -s,  --setup                  Set up whitelists for compatibility with new technology
  -t,  --technology PLATFORM    Name of technology used to generate data (10x, nadia, icell8)
  -R1, --read1 FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2 FILE             Read 2 FASTQ file to pass to cellranger
  -f,  --file NAME              Name of FASTQ files to pass to cellranger (prefix before R1 or R2)
  -b,  --barcodes FILE           custom iCELL8 barcode list in plain text
  -i,  --id ID                  A unique run id, used to name output folder
  -d,  --description TEXT       Sample description to embed in output files.
  -r,  --reference DIR          Path of directory containing 10x-compatible reference.
  -c,  --chemistry CHEM         Assay configuration, autodetection is not possible for converted files: 'SC3Pv2' (default), 'SC5P-PE', or 'SC5P-R2'
  -n,  --force-cells NUM        Force pipeline to use this number of cells, bypassing the cell detection algorithm.
  -j,  --jobmode MODE           Job manager to use. Valid options: 'local' (default), 'sge', 'lsf', or a .template file
  -p,  --pass                   Skips the FASTQ file conversion if converted files already exist
  -h,  --help                   Display this help and exit
  -v,  --version                Output version information and exit
       --verbose                Print additional outputs for debugging

For each fastq file, follow the naming convention below:
  <SampleName>_<SampleNumber>_<LaneNumber>_<ReadNumber>_001.fastq
  e.g. EXAMPLE_S1_L001_R1_001.fastq
       Example_S4_L002_R2_001.fastq.gz

Files will be renamed if they do not follow this format. File extension will be detected automatically.
'

if [[ -z $@ ]]; then
    echo "$help"
    exit 1
fi
##########



#####options#####
#set options
DIR=`which cellranger` #location of cellranger
VERSION=`cellranger count --version | head -n 2 | tail -n 1 | cut -d"(" -f2 | cut -d")" -f1` #get cellranger version
lockfile=${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.lock #path for .lock file
lastcallfile=${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.last_called #path for .last_called
lastcall=`cat $lastcallfile`
barcodefolder=${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes #folder with the barcodes
crIN=input4cellranger #name of the directory with all FASTQ files given to cellranger

#variable options
setup=false
convert=true
read1=()
read2=()
barcodes=""
SAMPLE=""
LANE=()
id=""
description=""
reference=""
ncells=""
chemistry=""
jobmode=""

next=false
for op in "$@"; do
    if $next; then
        next=false;
        continue;
    fi
    case "$op" in
        -v|--version)
            echo "$ver_info"
            exit 0
            ;;
        -h|--help)
            echo "$help"
            exit 0
            ;;
        -s|--setup)
            setup=true
            next=false
            shift
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
                echo "Error: File input missing --file or --read1"
                exit 1
            fi
            ;;
        -b|--barcodes)
            shift
            if [[ "$1" != "" ]]; then
                barcodes="${1/%\//}"
                next=true
                shift
            else
                echo "Error: value missing for --barcodes"
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
        -p|--pass)
            convert=false
            next=false
            shift
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



#####check if input maches expected inputs#####
if [[ $verbose == "true" ]]; then
    echo "checking options ..."
fi

#check if technology matches expected inputs
if [[ "$technology" != "10x" ]] && [[ "$technology" != "nadia" ]] && [[ "$technology" != "icell8" ]]; then
    echo "Error: option -t needs to be 10x, nadia, or icell8"
    exit 1
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
                rename "s/_$raadkey/_L001_$readkey/" $read
                #update file variable
                read=`echo $read | sed -e "s/_${readkey}/_L001_${raedkey}/g"`
                list[$j]=$read
            ;;
        esac
        case $read in
            #check if contains sample before lane
            *_S[123456789]_L0*)
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
                read=`echo $read | sed -e "s/_L0/_S${j}_L0/g"`
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
                read=`echo $read | sed -e "s/_${readkey}.*\./_${raedkey}_001\./g"`
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
    name=`basename $fq | cut -f1 -d'.' | grep -o "_" | wc -l | xargs`
    sn=`basename $fq | cut -f1-$(($name-3))  -d'_'`
    ln=`basename $fq | cut -f$(($name-1))  -d'_' | sed 's/L00//'`
    LANE+=($ln)
    if [[ $name < 4 ]]; then
        echo "Error: filename $fq is not following the naming convention. (e.g. EXAMPLE_S1_L001_R1_001.fastq)";
        exit 1
    elif [[ $fq != *'.fastq'* ]] && [[ $fq != *'.fq'* ]]; then
        echo "Error: $fq does not have a .fq or .fastq extention."
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

#checking the quality of custom barcode file
if [[ ! -z "$barcodes" ]]; then
    if [[ ! -f $barcodes ]]; then
        echo "Error: File selected for --barcode does not exist"
    elif [[ "$technology" != "icell8" ]]; then
        echo "Error: option --barcodes is only valid when --technology is icell8 yet its currently set as $technology."
        exit 1
    fi
fi
barcodes=`readlink -f $barcodes`

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

#check if chemistry matches expected input
if [[ -z "$chemistry" ]]; then
    chemistry="SC3Pv2"
elif [[ "$chemistry" != "SC3Pv2" ]] && [[ "$chemistry" != "SC5P-PE" ]] && [[ "$chemistry" != "SC5P-R2" ]]; then
    echo "Error: option --chemistry must be SC3Pv2, SC5P-PE , or SC5P-R2"
    exit 1
fi

#checking if jobmode matches expected input
if [[ -z "$jobmode" ]]; then
    jobmode="local"
elif [[ "$jobmode" != "local" ]] && [[ "$jobmode" != "sge" ]] && [[ "$jobmode" != "lsf" ]] && [[ "$jobmode" != *"template" ]]; then
    echo "Error: option --jobmode must be local, sge, lsf, or a .template file"
    exit 1
fi

#check if setup needs to be run before analysis (potentially overriding the user input)
if [[ -z $setup ]]; then
    setup=false
fi
if [[ $lastcall != $technology ]]; then
    setup=true
fi

#check if convertion is needs to be run before analysis (potentially overriding the user input)
if [[ ! -d $crIN ]] || [[ $lastcall != $technology ]]; then
    convert=true
fi

#check if ID is present
if [[ -z $id ]]; then
    if [[ ${#read1[@]} -ne 0 ]] || [[ ${#read2[@]} -ne 0 ]]; then
        echo "Error: option --id is required"
        exit 1
    fi
fi
crIN=${crIN}_${id}
if [[ ! -d ${crIN} ]]; then
    convert=true
fi
##########



#####check if UniverSC is running already#####
#set up .lock file
if [[ ! -f $lockfile ]]; then
    echo "creating .lock file"
    echo 0 > $lockfile
else
    #check if jobs running (check value in .lock file)
    echo "checking .lock file"
    lock=`cat $lockfile`
    
    if [[ $lock -eq 0 ]]; then
        echo " call accepted: no other cellranger jobs running"
        lock=1
        echo $lock > $lockfile
    else
        if [[ -f $lastcallfile ]]; then
            echo " total of $lock cellranger ${VERSION} jobs already running in ${DIR} with technology $lastcall"
            
            #check if the technology running is different from the current convert call
            if [[ $lastcall == "icell8_custom" ]]; then
                echo "Error: icell8 with custom barcode list is currently running"
                echo "other jobs cannot be run until the current job is complete"
                echo "remove $lockfile if $lastcall jobs have completed or aborted"
                exit 1
            elif [[ $lastcall == $technology ]]; then
                echo " call accepted: no conflict detected with other jobs currently running"
                #add current job to lock
                lock=$(($lock+1))
                echo $lock > $lockfile
                setup=false
            else
                echo "Error: conflict between technology selected for the new job ($technology) and other jobs currently running ($lastcall)"
                echo "barcode whitelist configured and locked for currently running technology: $lastcall"
                echo "remove $lockfile if $lastcall jobs have completed or aborted"
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
echo "SETUP: $setup"
echo "FORMAT: $technology"
if [[ $technology == "nadia" ]]; then
    echo "***Warning: whitelist is converted for compatibility with $technology, valid barcodes cannot be detected accurately with this technology***"
fi
if [[ -z $barcodes ]]; then
    echo "BARCODES: default"
else
    echo "BARCODES: (custom barcode file) $barcodes"
fi
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
echo "##########"
echo ""
##########



####setup whitelist#####
#run setup if called
if [[ $setup == "true" ]]; then
    echo "setup begin"
    echo "updating barcodes in $barcodefolder for cellranger version $VERSION installed in $DIR ..."
    
    if [[ ! -w $barcodefolder ]]; then
        echo "Error: launch_universc.sh can only be run cellranger installed locally"
        echo "Running cellranger installed at $DIR"
        echo "Install cellranger in a directory with write permissions such as /home/`whoami`/local"
        echo "cellranger must be exported to the PATH"
        echo "The following versions of cellranger are found:"
        echo " `whereis cellranger`"
        exit 1
    fi
    
    cd $barcodefolder
    
    #compressing v3 whilelist
    if [[ -f 3M-february-2018.txt ]]; then
        gzip -f 3M-february-2018.txt
    fi
    
    #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
    if [[ -f nadia_barcode.txt ]] || [[ -f iCell8_barcode.txt ]]; then
        echo " restoring 10x barcodes for version 2 kit ..."
        cp 737K-august-2016.txt.backup 737K-august-2016.txt
    fi
    echo " whitelist converted for 10x compatibility with version 2 kit."
    
    if [[ $technology == "10x" ]]; then
        #create version 3 files if version 3 whitelist available
        if [[ -f 3M-february-2018.txt.gz ]]; then
            #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
            if [[ -f nadia_barcode.txt ]] || [[ -f iCell8_barcode.txt ]]; then
                echo " restoring 10x barcodes for version 3 kit ..."
                cp 3M-february-2018.txt.backup.gz 3M-february-2018.txt.gz
            fi
            echo " whitelist converted for 10x compatibility with version 3 kit."
        fi
        
        #if cellranger version is 3 or greater, restore assert functions
        if [[ `printf '%s\n' '$VERSION 3.0.0' | sort -V | head -n 1` != $VERSION ]]; then
            #restore functions if other technology run previously
            if [[ $lastcall != $technology ]]; then
                sed -i "s/#if gem_group == prev_gem_group/if gem_group == prev_gem_group/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/#assert barcode_idx >= prev_barcode_idx/assert barcode_idx >= prev_barcode_idx/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/assert np.array_equal(in_mc.get_barcodes(), barcodes)/g" ${DIR}-cs/${VERSION}/lib/python/cellranger/molecule_counter.py
            fi
            echo " $DIR restored for $technology"
        else
            echo " $DIR ready for $technology"
        fi
    else
        #save original barcode file (if doesn't already exist)
        if [[ ! -f  737K-august-2016.txt.backup ]]; then
            echo " backing up whitelist of version 2 kit ..."
            cp 737K-august-2016.txt 737K-august-2016.txt.backup
        fi
        
        #create a new version 2 barcode file
        if [[ $technology == "nadia" ]]; then
            #create a Nadia barcode file
            if [[ ! -f nadia_barcode.txt ]]; then
                if [[ -f nadia_barcide.txt.gz ]]; then
                    gunzip -f nadia_barcode.txt.gz
                else
                    echo AAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > nadia_barcode.txt
                fi
            fi
            cat nadia_barcode.txt | sort | uniq > 737K-august-2016.txt
        elif [[ $technology == "icell8" ]]; then
            #create an iCell8 barcode file by copying from convert repo
            if [[ -f iCell8_barcode.txt.gz ]]; then
                rm iCell8_barcode.txt.gz
            fi
            if [[ -f iCell8_barcode.txt ]]; then
                rm iCell8_barcode.txt
            fi
            
            if [[ -z $barcodes ]]; then
                cat ${SCRIPT_DIR}/iCell8_barcode.txt >$barcodefolder/iCell8_barcode.txt
            else
                cat $barcodes >$barcodefolder/iCell8_barcode.txt
            fi
            sed -i 's/^/AAAAA/g' iCell8_barcode.txt #convert barcode whitelist to match converted barcodes
            cat iCell8_barcode.txt | sort | uniq > 737K-august-2016.txt
        else
            lock=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.lock`
            lock=$(($lock-1))
            echo $lock > $lockfile
            echo "Error: technology ($technology) is not supported"
            cd -
            exit 1
        fi
        echo " whitelist converted for $technology compatibility with version 2 kit."
        
        #create new version 3 files if version 3 whitelist available
        if [[ -f 3M-february-2018.txt.gz ]]; then
            #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
            if [[ -f nadia_barcode.txt ]] || [[ -f iCell8_barcode.txt ]]; then
                echo " restoring 10x barcodes for version 3 kit ..."
                cp 3M-february-2018.txt.backup.gz 3M-february-2018.txt.gz
            fi
            gunzip -k -f 3M-february-2018.txt.gz
            if [[ ! -f  3M-february-2018.txt.backup.gz ]]; then
                echo " backing up whitelist of version 3 kit ..."
                cp 3M-february-2018.txt 3M-february-2018.txt.backup
                gzip 3M-february-2018.txt.backup
            fi
            
            #combine 10x and new barcodes
            if [[ $technology == "nadia" ]]; then
                if [[ ! -f nadia_barcode.txt.gz ]]; then
                    gzip -f nadia_barcode.txt
                fi
                zcat nadia_barcode.txt.gz 3M-february-2018.txt.backup.gz | sort | uniq > 3M-february-2018.txt
            elif [[ $technology == "icell8" ]]; then
                if [[ -f iCell8_barcode.txt.gz ]]; then
                    rm iCell8_barcode.txt.gz
                fi
                if [[ -z $barcodes ]]; then
                    cat ${SCRIPT_DIR}/iCell8_barcode.txt > 3M-february-2018.txt
                else
                    cat ${barcodes} > 3M-february-2018.txt
                fi
            fi
            gzip -f 3M-february-2018.txt
            
            if [[ -d translation ]]; then
                rm translation/3M-february-2018.txt.gz
                ln -s ../3M-february-2018.txt.gz translation/3M-february-2018.txt.gz
            fi
        fi    
        echo " whitelist converted for $technology compatibility with version 3 kit."
        
        #if cellranger version is 3 or greater, restore assert functions
        if [[ `printf '%s\n' '$VERSION 3.0.0' | sort -V | head -n 1` != $VERSION ]]; then
            #disable functions if 10x technology run previously
            if [[ $lastcall == "10x" ]]; then
                sed -i "s/if gem_group == prev_gem_group/#if gem_group == prev_gem_group/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/assert barcode_idx >= prev_barcode_idx/#assert barcode_idx >= prev_barcode_idx/g" ${DIR}-cs/${VERSION}/mro/stages/counter/report_molecules/__init__.py
                sed -i "s/assert np.array_equal(in_mc.get_barcodes(), barcodes)/#assert np.array_equal(in_mc.get_barcodes(), barcodes)/g"  ${DIR}-cs/${VERSION}/lib/python/cellranger/molecule_counter.py
            fi
            echo " $DIR restored for $technology"
        else
            echo " $DIR ready for $technology"
        fi
    fi
    
    if [[ ! -z $barcodes ]]; then
        echo ${technology}_custom > $lastcallfile
    else
        echo $technology > $lastcallfile
    fi
    cd -
    echo "setup complete"
fi

if [[ ${#read1[@]} -eq 0 ]] && [[ ${#read2[@]} -eq 0 ]]; then
    lock=`cat $lockfile`
    lock=$(($lock-1))
    echo $lock > $lockfile
    echo " whitelist converted and no FASTQ files are selected. exiting launch_universc.sh"
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
    crR2s+=($to)
    
    echo " handling $fq ..."
    if [[ ! -f $to ]] || [[ $convert == "true" ]]; then
        cp -f $fq $to
    fi
done
##########



#####convert file format#####
echo "converting input files to confer cellranger format ..."
if [[ $convert = "false" ]]; then
    echo " input file format conversion skipped"
fi

if [[ "$technology" == "10x" ]]; then
    echo " 10x files accepted without conversion"
elif [[ $convert == "true" ]]; then
    echo " converting file from $technology format to 10x format ..."
    for convFile in "${convFiles[@]}"; do
        echo " handling $convFile ..."
        if [[ "$technology" == "nadia" ]]; then
            echo "  converting barcodes"
            sed -i '2~4s/^/AAAA/' $convFile #Add AAAA to every read
            echo "  converting quality scores"
            sed -i '4~4s/^/IIII/' $convFile #Add quality scores for added bases
            echo "  converting UMI"
            sed -i '2~4s/[NATCG][NATCG][NATCG][NATCG][NATCG][NATCG]$/AA/' $convFile #Replace last 6 bases with AA
            echo "  converting quality scores"
            sed -i '4~4s/......$/II/' $convFile #Replace quality scores for added bases
        elif [[ "$technology" == "icell8" ]]; then
            echo "  converting barcodes"
            sed -i '2~4s/^/AAAAA/' $convFile #Add AAAAA to every read
            echo "  converting quality scores"
            sed -i '4~4s/^/IIIII/' $convFile #Add quality scores for added bases
        fi
    done
fi
##########



#####run cellranger#####
echo "running cellranger ..."
echo ""
echo "#####cellranger#####"
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
echo "##########"
echo ""
##########



#####remove files if convert is not running elsewhere#####
echo "updating .lock file"

#remove currewnt job from counter (successfully completed)
lock=`cat ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes/.lock`
lock=$(($lock-1))

#check if jobs running
if [[ $lock -ge 1 ]]; then
    echo " total of $lock jobs for $lastcall technology are still run by cellranger ${VERSION} in ${DIR}"
else
    echo " no other jobs currently run by cellranger ${VERSION} in ${DIR}"
    echo " no conflicts: whitelist can now be changed for other technologies"
    rm -f $lockfile
fi
##########



#####printing out log#####
log="
#####Conversion tool log#####
$ver_info

Original barcode format: ${technology} (then converted to 10x)

Runtime: ${runtime}s
##########
"
echo "$log"
##########

exit 0
