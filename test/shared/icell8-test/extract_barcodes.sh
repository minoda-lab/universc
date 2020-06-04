#! /bin/bash
## usage: bash extract_barcodes.sh input.tsv output.txt
WellList=$1
BarcodeList=$2
tail -n $(($(wc -l $WellList | cut -d" " -f5)-1)) $WellList | cut -f6 > $BarcodeList
