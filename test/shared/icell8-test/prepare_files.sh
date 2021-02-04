#download files in conda environment (requires credentials)
conda activate ega  
pyega3 files EGAD00001003443
pyega3 -cf credential_file.json -d -t EGAD00001003443
pyega3 fetch EGAF00001698000
conda deactivate
mv EGAF00001698000/72618_KU812.fastq.gz 72618_KU812.fastq.gz

#export PATH
source $HOME/.bashrc  

#open compressed files
gunzip -fk 72618_KU812.fastq.gz
mv 72618_KU812.fastq 72618_KU812_R2_001.fastq

# create R1 from barcodes in headers
cat  72618_KU812_R2_001.fastq | sed -E '1~4s/ BC\:(.{11}) UMI\:(.{10})/ BC\:\1 UMI\:\2\n\1\2AAAAAAAAAAAAAAAAAAAAAAAAAAAAA/g' | sed '3~5d' > 72618_KU812_R1_001.fastq

#create bowtie index
bowtie2-build ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome
bowtie2 -x ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_R2_001.fastq -S 72618_KU812_R2_001.sam
samtools view -bS 72618_KU812_R2_001.sam > 72618_KU812_R2_001.bam

# extract region of genome
samtools sort -O BAM 72618_KU812_R2_001.bam > 72618_KU812_R2_001.sort.bam
samtools index 72618_KU812_R2_001.sort.bam
samtools view  72618_KU812_R2_001.sort.bam  21:9825832-48085036 > 72618_KU812_R2_001.chr21.bam
samtools view -O BAM  72618_KU812_R2_001.sort.bam  21:9825832-48085036 > 72618_KU812_R2_001.chr21.sort.bam
samtools sort -n 72618_KU812_R2_001.chr21.sort.bam -o 72618_KU812_R2_001.chr21.qsort.bam

# export mapped regions to fastq
bedtools bamtofastq -i 72618_KU812_R2_001.chr21.qsort.bam -fq 72618_KU812_R2_001.chr21_R2.fastq
fastq_pair 72618_KU812_R2_001.chr21_R2.fastq 72618_KU812_R1_001.fastq

# downsample to 250,000 reads per "lane"
seqtk sample -s999 72618_KU812_R1_001.fastq.paired.fq 250000 > 72618_KU812_L001_R1_001.fastq
seqtk sample -s999 72618_KU812_R2_001.chr21_R2.fastq.paired.fq 250000 > 72618_KU812_L001_R2_001.fastq
seqtk sample -s100 72618_KU812_R1_001.fastq.paired.fq 250000 > 72618_KU812_L002_R1_001.fastq
seqtk sample -s100 72618_KU812_R2_001.chr21_R2.fastq.paired.fq 250000 > 72618_KU812_L002_R2_001.fastq

# test mapping downsampled reads to genome
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L001_R2_001.fastq -S 72618_KU812_L001_R2_001.sam
samtools view -bS  72618_KU812_L001_R2_001.sam | samtools sort  >  72618_KU812_L001_R2_001.bam
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L002_R2_001.fastq -S 72618_KU812_L002_R2_001.sam
samtools view -bS  72618_KU812_L002_R2_001.sam | samtools sort  >  72618_KU812_L002_R2_001.bam

# sort and index
samtools sort -O BAM 72618_KU812_L001_R2_001.bam > 72618_KU812_L001_R2_001.sort.bam
samtools index 72618_KU812_L001_R2_001.sort.bam
samtools sort -O BAM 72618_KU812_L002_R2_001.bam > 72618_KU812_L002_R2_001.sort.bam
samtools index 72618_KU812_L002_R2_001.sort.bam

# check for variants in VCF format
bcftools mpileup -O b -o  72618_KU812_L001_R2_001.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L001_R2_001.sort.bam
bcftools mpileup -O b -o  72618_KU812_L002_R2_001.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L002_R2_001.sort.bam
bcftools call --ploidy 1 -m -v -o  72618_KU812_L001_R2_001.vcf 72618_KU812_L001_R2_001.bcf
bcftools call --ploidy 1 -m -v -o  72618_KU812_L002_R2_001.vcf 72618_KU812_L002_R2_001.bcf

# save headers
samtools view -H 72618_KU812_L001_R2_001.sam > 72618_KU812_L001_R2_001.header.sam
samtools view -H 72618_KU812_L002_R2_001.sam > 72618_KU812_L002_R2_001.header.sam

#call custom script to replace variants with reference sequence
Rscript update_sam_file_indels.R 72618_KU812_L001_R2_001.sam 72618_KU812_L001_R2_001.vcf
Rscript update_sam_file_indels.R 72618_KU812_L002_R2_001.sam 72618_KU812_L002_R2_001.vcf
# output: 72618_KU812_L001_R2_001.converted.sam

#remove trailing whitespace
sed 's/[ \t]*$//g' -i 72618_KU812_L001_R2_001.converted.sam 72618_KU812_L002_R2_001.converted.sam

# add headers
cat 72618_KU812_L001_R2_001.header.sam  72618_KU812_L001_R2_001.converted.sam >  72618_KU812_L001_R2_001.masked.sam
cat 72618_KU812_L002_R2_001.header.sam  72618_KU812_L002_R2_001.converted.sam >  72618_KU812_L002_R2_001.masked.sam

# export bam file
samtools view -bS 72618_KU812_L001_R2_001.masked.sam > 72618_KU812_L001_R2_001.masked.bam
samtools view -bS 72618_KU812_L002_R2_001.masked.sam > 72618_KU812_L002_R2_001.masked.bam

# sort and index
samtools sort -O BAM 72618_KU812_L001_R2_001.masked.bam > 72618_KU812_L001_R2_001.sort.masked.bam
samtools index 72618_KU812_L001_R2_001.sort.masked.bam
samtools sort -O BAM 72618_KU812_L002_R2_001.masked.bam > 72618_KU812_L002_R2_001.sort.masked.bam
samtools index 72618_KU812_L002_R2_001.sort.masked.bam

# export fastq
bedtools bamtofastq -i 72618_KU812_L001_R2_001.sort.masked.bam -fq 72618_KU812_L001_R2_001.masked.fastq
bedtools bamtofastq -i 72618_KU812_L002_R2_001.sort.masked.bam -fq 72618_KU812_L002_R2_001.masked.fastq

cp 72618_KU812_L001_R2_001.masked.fastq 72618_KU812_L001_R2_001.unmasked.fastq
cp 72618_KU812_L002_R2_001.masked.fastq 72618_KU812_L002_R2_001.unmasked.fastq

# match to read1 barcodes
cp 72618_KU812_L001_R1_001.fastq 72618_KU812_L001_R1_001.masked.fastq
fastq_pair 72618_KU812_L001_R2_001.masked.fastq 72618_KU812_L001_R1_001.masked.fastq
cp 72618_KU812_L002_R1_001.fastq 72618_KU812_L002_R1_001.masked.fastq
fastq_pair 72618_KU812_L002_R2_001.masked.fastq 72618_KU812_L002_R1_001.masked.fastq

cp 72618_KU812_L001_R2_001.masked.fastq.paired.fq 72618_KU812_L001_R2_001.unmasked.fastq.paired.fq
cp 72618_KU812_L002_R2_001.masked.fastq.paired.fq 72618_KU812_L002_R2_001.unmasked.fastq.paired.fq

# map corrected reads to genome
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L001_R2_001.masked.fastq.paired.fq -S 72618_KU812_L001_R2_001.check.masked.sam
samtools view -bS 72618_KU812_L001_R2_001.check.masked.sam > 72618_KU812_L001_R2_001.check.masked.bam
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L002_R2_001.masked.fastq.paired.fq -S 72618_KU812_L002_R2_001.check.masked.sam
samtools view -bS 72618_KU812_L002_R2_001.check.masked.sam > 72618_KU812_L002_R2_001.check.masked.bam

# sort and index
samtools sort -O BAM 72618_KU812_L001_R2_001.check.masked.bam > 72618_KU812_L001_R2_001.sort.check.masked.bam
samtools index 72618_KU812_L001_R2_001.sort.check.masked.bam
samtools sort -O BAM 72618_KU812_L002_R2_001.check.masked.bam > 72618_KU812_L002_R2_001.sort.check.masked.bam
samtools index 72618_KU812_L002_R2_001.sort.check.masked.bam

# check for variants in VCF format
bcftools mpileup -O b -o  72618_KU812_L001_R2_001.masked.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L001_R2_001.sort.check.masked.bam
bcftools mpileup -O b -o  72618_KU812_L002_R2_001.masked.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L002_R2_001.sort.check.masked.bam
bcftools call --ploidy 1 -m -v -o  72618_KU812_L001_R2_001.masked.vcf 72618_KU812_L001_R2_001.masked.bcf
bcftools call --ploidy 1 -m -v -o  72618_KU812_L002_R2_001.masked.vcf 72618_KU812_L002_R2_001.masked.bcf

cp  72618_KU812_L001_R2_001.masked.vcf  72618_KU812_L001_R2_001.vcf
cp  72618_KU812_L002_R2_001.masked.vcf  72618_KU812_L002_R2_001.vcf
cp  72618_KU812_L001_R2_001.masked.sam  72618_KU812_L001_R2_001.sam
cp  72618_KU812_L002_R2_001.masked.sam  72618_KU812_L002_R2_001.sam

#call custom script to replace variants with reference sequence
Rscript update_sam_file_snps.R 72618_KU812_L001_R2_001.sam 72618_KU812_L001_R2_001.vcf
Rscript update_sam_file_snps.R 72618_KU812_L002_R2_001.sam 72618_KU812_L002_R2_001.vcf
# output: 72618_KU812_L001_R2_001.converted.sam

#remove trailing whitespace
sed 's/[ \t]*$//g' -i 72618_KU812_L001_R2_001.converted.sam 72618_KU812_L002_R2_001.converted.sam

# add headers
cat 72618_KU812_L001_R2_001.header.sam  72618_KU812_L001_R2_001.converted.sam >  72618_KU812_L001_R2_001.masked.sam
cat 72618_KU812_L002_R2_001.header.sam  72618_KU812_L002_R2_001.converted.sam >  72618_KU812_L002_R2_001.masked.sam

# export bam file
samtools view -bS 72618_KU812_L001_R2_001.masked.sam > 72618_KU812_L001_R2_001.masked.bam
samtools view -bS 72618_KU812_L002_R2_001.masked.sam > 72618_KU812_L002_R2_001.masked.bam

# sort and index
samtools sort -O BAM 72618_KU812_L001_R2_001.masked.bam > 72618_KU812_L001_R2_001.sort.masked.bam
samtools index 72618_KU812_L001_R2_001.sort.masked.bam
samtools sort -O BAM 72618_KU812_L002_R2_001.masked.bam > 72618_KU812_L002_R2_001.sort.masked.bam
samtools index 72618_KU812_L002_R2_001.sort.masked.bam

# export fastq
bedtools bamtofastq -i 72618_KU812_L001_R2_001.sort.masked.bam -fq 72618_KU812_L001_R2_001.masked.fastq
bedtools bamtofastq -i 72618_KU812_L002_R2_001.sort.masked.bam -fq 72618_KU812_L002_R2_001.masked.fastq

# match to read1 barcodes
cp 72618_KU812_L001_R1_001.fastq 72618_KU812_L001_R1_001.masked.fastq
fastq_pair 72618_KU812_L001_R2_001.masked.fastq 72618_KU812_L001_R1_001.masked.fastq
cp 72618_KU812_L002_R1_001.fastq 72618_KU812_L002_R1_001.masked.fastq
fastq_pair 72618_KU812_L002_R2_001.masked.fastq 72618_KU812_L002_R1_001.masked.fastq

# map corrected reads to genome
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L001_R2_001.masked.fastq.paired.fq -S 72618_KU812_L001_R2_001.check.masked.sam
samtools view -bS 72618_KU812_L001_R2_001.check.masked.sam > 72618_KU812_L001_R2_001.check.masked.bam
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U 72618_KU812_L002_R2_001.masked.fastq.paired.fq -S 72618_KU812_L002_R2_001.check.masked.sam
samtools view -bS 72618_KU812_L002_R2_001.check.masked.sam > 72618_KU812_L002_R2_001.check.masked.bam

# sort and index
samtools sort -O BAM 72618_KU812_L001_R2_001.check.masked.bam > 72618_KU812_L001_R2_001.sort.check.masked.bam
samtools index 72618_KU812_L001_R2_001.sort.check.masked.bam
samtools sort -O BAM 72618_KU812_L002_R2_001.check.masked.bam > 72618_KU812_L002_R2_001.sort.check.masked.bam
samtools index 72618_KU812_L002_R2_001.sort.check.masked.bam

# check for variants in VCF format
bcftools mpileup -O b -o  72618_KU812_L001_R2_001.masked.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L001_R2_001.sort.check.masked.bam
bcftools mpileup -O b -o  72618_KU812_L002_R2_001.masked.bcf -f ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa 72618_KU812_L002_R2_001.sort.check.masked.bam
bcftools call --ploidy 1 -m -v -o  72618_KU812_L001_R2_001.masked.vcf 72618_KU812_L001_R2_001.masked.bcf
bcftools call --ploidy 1 -m -v -o  72618_KU812_L002_R2_001.masked.vcf 72618_KU812_L002_R2_001.masked.bcf

rename 's/sort.masked/masked/g' *bam
rename -f 's/.masked//g' *.masked*
rename -f 's/\.paired\.fq//g' *.paired.fq
