#! /bin/bash

cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
cellrangerpath=`which cellranger`
cellrangerdir=$(dirname $(which cellranger))

#echo "check for SA files"
echo $cellrangerdir $cellrangerpath $cellrangerversion
if [[ ! -f 3.0.0/star/SA ]] && [[ -f ${cellrangerpathdir}/cellranger-tiny-ref/3.0.0/star/SA ]]; then
    rsync ${cellrangerdir}/cellranger-tiny-ref/3.0.0/star/SA 3.0.0/star/SA
fi
if [[ ! -f 1.2.0/star/SA ]] && [[ -f ${cellrangerdir}/cellranger-tiny-ref/1.2.0/star/SA ]]; then
    rsync ${cellrangerdir}/cellranger-tiny-ref/1.2.0/star/SA 1.2.0/star/SA
fi
touch 3.0.0/star/SA
touch 1.2.0/star/SA
