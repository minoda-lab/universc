#find barcodes directory of local cellranger install
#cd /home/tom/local/bin/cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes
#cd /home/tom/local/bin/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes
DIR=`which /home/tom/local/bin/cellranger-2.1.0/cellranger`
VERSION=`cellranger count --version | head -n 2 | tail -n 1 | cut -d"(" -f2 | cut -d")" -f1`
cd ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes

echo "update barcodes in ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes \n for cellranger version $VERSION installed in $DIR"

#restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
    then
    echo "restore 10x barcodes"
    cp 737K-august-2016.txt.backup 737K-august-2016.txt
    fi

#create a file with every possible barcode (permutation)
if [ ! -f nadia_barcode.txt ]
    then
    echo AAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > nadia_barcode.txt
    echo "expected barcodes generated for Nadia"
    fi 

#save original barcode file (if doesn't already exist)
if [ ! -f  737K-august-2016.txt.backup ]
    then
    echo backup of version 2 whitelist
    cp 737K-august-2016.txt 737K-august-2016.txt.backup
    fi

#combine 10x and Nadia barcodes
cat nadia_barcode.txt 737K-august-2016.txt.backup > 737K-august-2016.txt
echo "whitelist converted for Nadia compatibility with version 2 kit"

#create version 3 files if version 3 whitelist available
if [ -f 3M-february-2018.txt.gz ]
    then
    #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
    if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
        then
        echo "restore 10x barcodes"
        cp 3M-february-2018.txt.gz.backup 3M-february-2018.txt.gz
        fi
    gunzip -k  3M-february-2018.txt.gz
    if [ ! -f  3M-february-2018.txt.backup.gz ]
        then
        echo "backup of version 3 whitelist"
        cp 3M-february-2018.txt 3M-february-2018.txt.backup
        gzip 3M-february-2018.txt.backup
        fi
    #combine 10x and Nadia barcodes
    gzip -k nadia_barcode.txt
    zcat nadia_barcode.txt 3M-february-2018.txt.backup > 3M-february-2018.txt
    gzip -f 3M-february-2018.txt
    echo "whitelist converted for Nadia compatibility with version 3 kit"
fi

