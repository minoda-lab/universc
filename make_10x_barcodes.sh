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
    echo "restore 10x barcodes
    cp 737K-august-2016.txt.backup 737K-august-2016.txt
    fi

echo "whitelist converted for 10x compatibility with version 2 kit"

#create version 3 files if version 3 whitelist available
if [ -f 3M-february-2018.txt.gz ]
    then
    #restore 10x barcodes if scripts has already been run (allows changing Nadia to iCELL8)
    if [ -f nadia_barcode.txt -o -f  iCELL8_barcode.txt ]
        then
        echo "restore 10x barcodes
        cp 3M-february-2018.txt.gz.backup 3M-february-2018.txt.gz
        fi
    echo "whitelist converted for 10x compatibility with version 3 kit"
fi

