#!/bin/bash


######convert version#####
crconverterversion="3.1.1.90001"
##########



#####cellranger version#####
cellrangerpath=`which cellranger` #location of Cell Ranger
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
of Cell Ranger and may not produce the same results as
software produced by 10x Genomics. We recommend the
supported software be used for other purposes.

The code in this repository is licensed under the MIT License.
'
##########




#####print usage#####
if [[ -z $1 ]]; then
    echo "$statement"
    exit 1
fi
##########




#####declare variables#####
name=$1
type=$2
inpath=""
output=""
pipestance=""
metrics=""
analysis=""
matrix=""
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

if [[ $1 == "-h" ]] || [[ $2 == "-h" ]] || [[ $1 == "--help" ]] || [[ $2 == "--help" ]]; then
    end=true
    usage=true
    help=true
fi

if [[ $1 == "-v" ]] || [[ $2 == "-v" ]] || [[ $1 == "--version" ]] || [[ $2 == "--version" ]]; then
    end=true
    version=true
fi
shift
shift
for op in "${@}"; do
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
                shift
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
          --matrix)
            shift
            if [[ "$1" != "" ]]; then
                matrix="${1/%\//}"
                next=true
                shift
            else
                matrix=""
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
                description=""
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
            echo "Error: Invalid option: $op"
            exit 1
            ;;
    esac
done
##########

if [[ $end  ]]; then
    if [[ ! -z $help ]]; then
        echo "$statement"
        echo "$options"
        exit 0
    fi
fi

if [[ $end ]]; then
    if [[ ! -z $version ]]; then
        echo "$crconverterversion"
        exit 0
    fi
fi

##########

if [[ -z $output ]]; then
    echo "An output file must be specified"
    exit 1
fi

if [[ ! -f $output ]]; then
    echo "#!/bin/bash" > $output
    echo "echo \"Please use the supported release of Cell Ranger if you wish to use .cloupe files\"" >> $output
    chmod 755 $output
else
    echo ".cloupe file already exists in $output"
fi

echo "Wrote .cloupe file to $output"
echo "Warning .cloupe file may be dysfunctional"
exit 0

