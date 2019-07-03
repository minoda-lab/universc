#find barcodes directory of local cellranger install
#cd /home/tom/local/bin/cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes
#cd /home/tom/local/bin/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes
DIR=`which /home/tom/local/bin/cellranger-2.1.0/cellranger`
VERSION=`cellranger count --version | head -n 2 | tail -n 1 | cut -d"(" -f2 | cut -d")" -f1`
cd ${DIR}-cs/${VERSION}/lib/python/cellranger/barcodes

#create a file with every possible barcode (permutation)
echo AAAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > iCELL8_barcode.txt

#save original barcode file (if doesn't already exist)
if [ ! -f  737K-august-2016.txt.backup ] then
    cp 737K-august-2016.txt 737K-august-2016.txt.backup
    fi

#combine 10x and Nadia barcodes
cat iCELL8_barcode.txt 737K-august-2016.txt.backup > 737K-august-2016.txt
