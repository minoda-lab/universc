#!/bin/bash


######convert version#####
crconverterversion="3.1.1.90001"
##########



#####cellranger version#####
cellrangerpath=`which cellranger` #location of cellranger
if [[ -z $cellrangerpath ]]; then
    echo "cellranger command is not found."
    exit 1
fi
cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
##########



#####declare variables#####
SOURCE="${BASH_SOURCE[0]}"
##########



#####usage statement#####
statement='crconverter

 Usage:

        crconverter <name> <type> [options]
    crconverter -h | --help | --version
'
options='

Options:
        --inpath=PATH          Path to a completed pipestance. crconvert will attempt find inputs automatically by filename
                               Explicit setting of input using other flags will override inputs found with --inpath
    --output=PATH      
    --pipestance=PATH 
    --metrics=PATH     
    --analysis=PATH    
    --matrix=PATH      
    --aggregation=PATH 
        --gemgroups=PATH   
        --contiginfo=PATH
        --peaks=PATH
        --fragmentsindex=PATH
        --geneannotations=PATH
        --geneannotationtypes=TYPES
    --description=DESC 
    -h --help          
    --version          

This is an open-source implementation of crconverter.
This is only intended for testing open-source implementations
of cellranger and may not produce the same results as
software produced by 10x Genomics. We recommend the
supported software be used for other purposes.

The code in this repository is licensed under the MIT License.
'

##########

#variable optionsname=false
name=$1
type=$2
inpath=""
output=""
pipestance=""
metrics=""
analysis=""
aggregation=""
gemgroups=""
contiginfo=""
peaks=""
fragmentsindex=""
geneannotations=""
geneannotationtypes=""
description=""

next=false
end=false

# crconverter tiny SC_RNA_COUNTER_CS --matrix /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/_BASIC_SC_RNA_COUNTER/FILTER_BARCODES/fork0/join-u00fa679091/files/filtered_matrices_h5.h5 \
# --analysis /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/SC_RNA_ANALYZER/SUMMARIZE_ANALYSIS/fork0/join-u00fa67911f/files/analysis/analysis.h5 \ 
# --output /root/tiny/SC_RNA_COUNTER_CS/CLOUPE_PREPROCESS/fork0/chnk0-u3dd5684351/files/output_for_cloupe.cloupe \ 
# --description CellRangerTest \ 
# --metrics /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/SUMMARIZE_REPORTS/fork0/join-u00fa679151/files/metrics_summary_json.json \
# --gemgroups /root/tiny/SC_RNA_COUNTER_CS/CLOUPE_PREPROCESS/fork0/chnk0-u3dd5684351/files/gem_group_index_json.json




for op in "${@:3}"; do
    if $next; then
        next=false;
        continue;
    fi
    case "$op" in
          --inpath)
            shift
            if [[ "$1" != "" ]]; then
                inpath="${1/%\//}"
                next=true
                shift
            else
                echo "Explicit setting of input using other flags will override inputs found with --inpath"
            fi
            ;;
           --output)
            shift
            if [[ "$1" != "" ]]; then
                output="${1/%\//}"
                next=true
                shift
            else
                output=""
            fi
            ;;
           --pipestance)
            shift
            if [[ "$1" != "" ]]; then
                pipestance="${1/%\//}"
                next=true
                shift
            else
                pipestance=""
            fi
            ;;
           --metrics)
            shift
            if [[ "$1" != "" ]]; then
                metrics="${1/%\//}"
                next=true
                shift
            else
                metrics=""
            fi
            ;;
           --analysis)
            shift
            if [[ "$1" != "" ]]; then
                analysis="${1/%\//}"
                next=true
                shift
            else
                analysis=""
            fi
            ;;
           --aggregation)
            shift
            if [[ "$1" != "" ]]; then
                aggregation="${1/%\//}"
                next=true
                shift
            else
                aggregation=""
            fi
            ;;
           --gemgroups)
            shift
            if [[ "$1" != "" ]]; then
                gemgroups="${1/%\//}"
                next=true
                shift
            else
                gemgroups=""
            fi
            ;;
           --contiginfo)
            shift
            if [[ "$1" != "" ]]; then
                contiginfo="${1/%\//}"
                next=true
                shift
            else
                contiginfo=""
            fi
            ;;
           --peaks)
            shift
            if [[ "$1" != "" ]]; then
                peaks="${1/%\//}"
                next=true
                shift
            else
                peaks=""
            fi
            ;;
           --fragmentsindex)
            shift
            if [[ "$1" != "" ]]; then
                fragmentsindex="${1/%\//}"
                next=true
                shift
            else
                fragmentsindex=""
            fi
            ;;
           --geneannotations)
            shift
            if [[ "$1" != "" ]]; then
                geneannotations="${1/%\//}"
                next=true
                shift
            else
                geneannotations=""
            fi
            ;;
           --geneannotationtypes)
            shift
            if [[ "$1" != "" ]]; then
                geneannotationtypes="${1/%\//}"
                next=true
                shift
            else
                geneannotationtypes=""
            fi
            ;;
           --description)
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
        -v|--version)
            end=true
            version=true
            echo "$crconverterversion"
            exit 0
            ;;
        -h|--help)
            end=true
            usage=true
            help=true
            echo "$statement"
            echo "$options"
            exit 0
            ;;
       --verbose)
            echo "debugging mode activated"
            verbose=true
            next=false
            shift
            ;;
        -*)
            echo "$statement"
            exit 1
            ;;
    esac
done
##########

if [[ -z $1 ]]; then
    echo "$statement"
    exit 1
fi

if [[ -z $help ]]; then
    echo "$statement"
    echo "$options"
    exit 0
fi

if [[ -z $version ]]; then
    echo "$crconverterversion"
    exit 0
fi

##########

echo "Cloupe File not generated"
exit 0

#Wrote .cloupe file to /home/tom/datasets/20190717_Plant_Protoplast_fiveprime_HiSeq/test_AtRTD2_P09_PE/SC_RNA_COUNTER_CS/CLOUPE_PREPROCESS/fork0/chnk0-u4b566787c3/files/output_for_cloupe.cloupe

# crconverter tiny SC_RNA_COUNTER_CS --matrix /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/_BASIC_SC_RNA_COUNTER/FILTER_BARCODES/fork0/join-u00fa679091/files/filtered_matrices_h5.h5 \
# --analysis /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/SC_RNA_ANALYZER/SUMMARIZE_ANALYSIS/fork0/join-u00fa67911f/files/analysis/analysis.h5 \ 
# --output /root/tiny/SC_RNA_COUNTER_CS/CLOUPE_PREPROCESS/fork0/chnk0-u3dd5684351/files/output_for_cloupe.cloupe \ 
# --description CellRangerTest \ 
# --metrics /root/tiny/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/SUMMARIZE_REPORTS/fork0/join-u00fa679151/files/metrics_summary_json.json \
# --gemgroups /root/tiny/SC_RNA_COUNTER_CS/CLOUPE_PREPROCESS/fork0/chnk0-u3dd5684351/files/gem_group_index_json.json


