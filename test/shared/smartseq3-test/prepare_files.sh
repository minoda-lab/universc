wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.R1.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.R2.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.I1.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.I2.fastq.gz
rename "s/Smartseq3.Fibroblasts.GelCut./Smartseq3_Fibroblasts_GelCut_/g" FSmartseq3.Fibroblasts.GelCut.*fastq.gz

#create bowtie index
bowtie2-build ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome
bowtie2 -x ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U Smartseq3_Fibroblasts_GelCut_R2_001.fastq -S Smartseq3_Fibroblasts_GelCut_R2_001.sam
samtools view -bS Smartseq3_Fibroblasts_GelCut_R2_001.sam > Smartseq3_Fibroblasts_GelCut_R2_001.bam

# extract region of genome
samtools sort -O BAM Smartseq3_Fibroblasts_GelCut_R2_001.bam > Smartseq3_Fibroblasts_GelCut_R2_001.sort.bam
samtools index Smartseq3_Fibroblasts_GelCut_R2_001.sort.bam
samtools view  Smartseq3_Fibroblasts_GelCut_R2_001.sort.bam  21:9825832-48085036 > Smartseq3_Fibroblasts_GelCut_R2_001.chr21.bam
samtools view -O BAM  Smartseq3_Fibroblasts_GelCut_R2_001.sort.bam  21:9825832-48085036 > Smartseq3_Fibroblasts_GelCut_R2_001.chr21.sort.bam
samtools sort -n Smartseq3_Fibroblasts_GelCut_R2_001.chr21.sort.bam -o Smartseq3_Fibroblasts_GelCut_R2_001.chr21.qsort.bam

# export mapped regions to fastq
bedtools bamtofastq -i Smartseq3_Fibroblasts_GelCut_R2_001.chr21.qsort.bam -fq Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq
fastq_pair Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq Smartseq3_Fibroblasts_GelCut_R1_001.fastq
fastq_pair Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq Smartseq3_Fibroblasts_GelCut_I1_001.fastq
fastq_pair Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq Smartseq3_Fibroblasts_GelCut_I2_001.fastq

# downsample to 250,000 reads per "lane"
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R1_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_R1_001.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_R2_001.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I1_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_I1_001.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I2_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_I2_001.fastq
seqtk sample -s100 Smartseq3_Fibroblasts_GelCut_R1_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L002_R1_001.fastq
seqtk sample -s100 Smartseq3_Fibroblasts_GelCut_R2_001.chr21_R2.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L002_R2_001.fastq
seqtk sample -s100 Smartseq3_Fibroblasts_GelCut_I1_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_I1_001.fastq
seqtk sample -s100 Smartseq3_Fibroblasts_GelCut_I2_001.fastq.paired.fq 250000 > Smartseq3_Fibroblasts_GelCut_L001_I2_001.fastq

#subsample smaller files for testing
rename "s/Smartseq3.Fibroblasts.GelCut./Smartseq3_Fibroblasts_GelCut_/g" Smartseq3.Fibroblasts.GelCut.*fastq*
# subsample to 500,000 reads
zcat Smartseq3_Fibroblasts_GelCut_R1.fastq.gz | head -n 2000000 > Smartseq3_Fibroblasts_GelCut_small_R1.fastq
zcat Smartseq3_Fibroblasts_GelCut_R2.fastq.gz | head -n 2000000 > Smartseq3_Fibroblasts_GelCut_small_R2.fastq
zcat Smartseq3_Fibroblasts_GelCut_I1.fastq.gz | head -n 2000000 > Smartseq3_Fibroblasts_GelCut_small_I1.fastq
zcat Smartseq3_Fibroblasts_GelCut_I2.fastq.gz | head -n 2000000 > Smartseq3_Fibroblasts_GelCut_small_I2.fastq
# subsample to 2.5 million reads
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R1.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_med_R1.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R2.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_med_R2.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I1.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_med_I1.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I2.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_med_I2.fastq
# subsample to 25 million reads
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R1.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_large_R1.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_R2.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_large_R2.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I1.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_large_I1.fastq
seqtk sample -s999 Smartseq3_Fibroblasts_GelCut_I2.fastq.gz 2500000 > Smartseq3_Fibroblasts_GelCut_large_I2.fastq

