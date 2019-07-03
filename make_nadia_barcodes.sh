#find barcodes directory of local cellranger install
cd /home/tom/local/bin/cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes

#create a file with every possible barcode (permutation)
echo AAAA{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G}{A,T,C,G} | sed 's/ /\n/g' > nadia_barcode.txt

#save original barcode file
cp 737K-august-2016.txt 737K-august-2016.txt.backup

#combine 10x and Nadia barcodes
cat nadia_barcodes.txt 737K-august-2016.txt.backup > 737K-august-2016.txt
