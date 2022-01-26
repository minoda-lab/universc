wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.R1.fastq.gz .
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.R2.fastq.gz .
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.I1.fastq.gz .
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.I2.fastq.gz . 
  
rename "s/.fastq/_001.fastq/g" Smartseq3_diySpike_S1*gz
rename "s/.diySpike./_diySpike_S1_L001_/g" *gz

#create bowtie index
bowtie2-build ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome
bowtie2 -x ../../test/cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U Smartseq3_diySpike_S1_L001_R2_001.fastq -S Smartseq3_diySpike_S1_L001_R2_001.sam
samtools view -bS Smartseq3_diySpike_S1_L001_R2_001.sam > Smartseq3_diySpike_S1_L001_R2_001.bam

# extract region of genome
samtools sort -O BAM Smartseq3_diySpike_S1_L001_R2_001.bam > Smartseq3_diySpike_S1_L001_R2_001.sort.bam
samtools index Smartseq3_diySpike_S1_L001_R2_001.sort.bam
samtools view  Smartseq3_diySpike_S1_L001_R2_001.sort.bam  21:9825832-48085036 > Smartseq3_diySpike_S1_L001_R2_001.chr21.bam
samtools view -O BAM  Smartseq3_diySpike_S1_L001_R2_001.sort.bam  21:9825832-48085036 > Smartseq3_diySpike_S1_L001_R2_001.chr21.sort.bam
samtools sort -n Smartseq3_diySpike_S1_L001_R2_001.chr21.sort.bam -o Smartseq3_diySpike_S1_L001_R2_001.chr21.qsort.bam

# export mapped regions to fastq
bedtools bamtofastq -i Smartseq3_diySpike_S1_L001_R2_001.chr21.qsort.bam -fq Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq
fastq_pair Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq Smartseq3_diySpike_S1_L001_R1_001.fastq
fastq_pair Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq Smartseq3_diySpike_S1_L001_I1_001.fastq
fastq_pair Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq Smartseq3_diySpike_S1_L001_I2_001.fastq

# downsample to 250,000 reads per "lane"
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R1_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L001_R1_001.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L001_R2_001.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I1_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L001_I1_001.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I2_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L001_I2_001.fastq
seqtk sample -s100 Smartseq3_diySpike_S1_L001_R1_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L002_R1_001.fastq
seqtk sample -s100 Smartseq3_diySpike_S1_L001_R2_001.chr21_R2.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L002_R2_001.fastq
seqtk sample -s100 Smartseq3_diySpike_S1_L001_I1_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L002_I1_001.fastq
seqtk sample -s100 Smartseq3_diySpike_S1_L001_I2_001.fastq.paired.fq 250000 > test-smartseq-hek293t_S2_L002_I2_001.fastq

#subsample smaller files for testing
rename "s/Smartseq3_diySpike_S1_L001_/Smartseq3_diySpike_S1_L001_/g" Smartseq3_diySpike_S1_L001_*fastq*
# subsample to 500,000 reads
zcat Smartseq3_diySpike_S1_L001_R1_001.fastq.gz | head -n 2000000 > Smartseq3_diySpike_S1_L001_small_R1.fastq
zcat Smartseq3_diySpike_S1_L001_R2_001.fastq.gz | head -n 2000000 > Smartseq3_diySpike_S1_L001_small_R2.fastq
zcat Smartseq3_diySpike_S1_L001_I1_001.fastq.gz | head -n 2000000 > Smartseq3_diySpike_S1_L001_small_I1.fastq
zcat Smartseq3_diySpike_S1_L001_I2_001.fastq.gz | head -n 2000000 > Smartseq3_diySpike_S1_L001_small_I2.fastq
# subsample to 2.5 million reads
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R1_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_med_R1.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R2_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_med_R2.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I1_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_med_I1.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I2_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_med_I2.fastq
# subsample to 25 million reads
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R1_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_large_R1.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_R2_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_large_R2.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I1_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_large_I1.fastq
seqtk sample -s999 Smartseq3_diySpike_S1_L001_I2_001.fastq.gz 2500000 > Smartseq3_diySpike_S1_L001_large_I2.fastq

