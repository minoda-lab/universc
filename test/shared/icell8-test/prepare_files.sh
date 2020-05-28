wget https://www.ncbi.nlm.nih.gov/geo/download/\?acc\=${1}\&format\=file\&file\=${1}%5FPure%5FHumanMouse%2Ebam
mv index.html\?acc=${1}\&format=file\&file=${1}%5FPure%5FHumanMouse%2Ebam ${1}.bam
samtools sort -n ${1}.bam > ${1}.qsort
samtools view  ${1}.qsort  21:9825832-48085036 > ${1}.qsort2
samtools sort -O BAM ${1}.bam > ${1}.sort.bam
samtools index ${1}.sort.bam
samtools view  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.bam
samtools view -O BAM  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.sort.bam
samtools sort -n ${1}.chr21.sort.bam -o ${1}.chr21.qsort.bam
bedtools bamtofastq -i ${1}.chr21.qsort.bam -fq ${1}_chr21_R1.fastq
mv ${1}_chr21_R1.fastq ${1}_chr21_R2.fastq
fastq-dump -F --split-files SRR1873277
fastq_pair ${1}_chr21_R2.fastq SRR1873277_1.fastq
head -n 117060 SRR1873277_1.fastq.paired.fq 117060 > SRR1873277_1.fastq.paired.fq
head -n 117060 ${1}_chr21_R2.fastq.paired.fq > ${1}_chr21_R2.fastq.paired.fq
cp SRR1873277_1.fastq.paired.fq  ${1}_chr21_R2.fastq.paired.fq ~/repos/universc/test/shared/dropseq-test
cp SRR1873277_1.fastq.paired.fq  ${1}_chr21_R2.fastq.paired.fq ~/repos/universc/test/shared/dropseq-test
mv SRR1873277_1.fastq.paired.fq SRR1873277_R1.fastq
mv ${1}_chr21_R2.fastq.paired.fq  SRR1873277_R2.fastq
mv ${1}_chr21_R2.fastq.paired.fq  SRR1873277_R2.fastq

bowtie2-build ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome.fa ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome

bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U iCELL8_01_S1_L001_R2_001.fastq | samtools sort > iCELL8_01_S1_L001_R2_001.bam
bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U iCELL8_01_S1_L002_R2_001.fastq | samtools sort > iCELL8_01_S1_L002_R2_001.bam

1=iCELL8_01_S1_L001_R2_001
samtools view -bS ${1}.sam > ${1}.bam

samtools sort -n ${1}.bam > ${1}.qsort
samtools view  ${1}.qsort  21:9825832-48085036 > ${1}.qsort2
samtools sort -O BAM ${1}.bam > ${1}.sort.bam
samtools index ${1}.sort.bam
samtools view  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.bam
samtools view -O BAM  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.sort.bam
samtools sort -n ${1}.chr21.sort.bam -o ${1}.chr21.qsort.bam
bedtools bamtofastq -i ${1}.chr21.qsort.bam -fq ${1}_chr21_R1.fastq
mv ${1}_chr21_R1.fastq ${1}_chr21_R2.fastq

1=iCELL8_01_S1_L002_R2_001
samtools sort -n ${1}.bam > ${1}.qsort
samtools view  ${1}.qsort  21:9825832-48085036 > ${1}.qsort2
samtools sort -O BAM ${1}.bam > ${1}.sort.bam
samtools index ${1}.sort.bam
samtools view  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.bam
samtools view -O BAM  ${1}.sort.bam  21:9825832-48085036 > ${1}.chr21.sort.bam
samtools sort -n ${1}.chr21.sort.bam -o ${1}.chr21.qsort.bam
bedtools bamtofastq -i ${1}.chr21.qsort.bam -fq ${1}_chr21_R1.fastq
mv ${1}_chr21_R1.fastq ${1}_chr21_R2.fastq



samtools index  iCELL8_01_S1_L001_R2_001.bam   
samtools index  iCELL8_01_S1_L002_R2_001.bam  

samtools sort -O BAM iCELL8_01_S1_L001_R2_001.bam > iCELL8_01_S1_L001_R2_001.sort.bam
samtools sort -O BAM iCELL8_01_S1_L002_R2_001.bam > iCELL8_01_S1_L002_R2_001.sort.bam
samtools index  iCELL8_01_S1_L001_R2_001.sort.bam   
samtools index  iCELL8_01_S1_L002_R2_001.sort.bam  
samtools view iCELL8_01_S1_L001_R2_001.sort.bam  21:9825832-48085036 > iCELL8_01_S1_L001_R2_001.chr21.bam    
samtools view iCELL8_01_S1_L002_R2_001.sort.bam  21:9825832-48085036 > iCELL8_01_S1_L002_R2_001.chr21.bam 

bowtie2 -x ../../cellranger_reference/cellranger-tiny-ref/3.0.0/fasta/genome -U /home/tom/datasets/20190408_iCELL8_MCF10A/iCELL8_01_S1_L001_R2_001.fastq -S iCELL8_01_S1_L001_R2_001.sam
samtools view -bS iCELL8_01_S1_L001_R2_001.sam > iCELL8_01_S1_L001_R2_001.bam
samtools sort -n iCELL8_01_S1_L001_R2_001.bam > iCELL8_01_S1_L001_R2_001.qsort

fastq_pair /home/tom/datasets/20190408_iCELL8_MCF10A/iCELL8_01_S1_L001_R1_001.fastq iCELL8_01_S1_L001_R2_001_2_chr21_R2.fastq
fastq_pair iCELL8_01_S1_L001_R2_001_2_chr21_R2.fastq /home/tom/datasets/20190408_iCELL8_MCF10A/iCELL8_01_S1_L001_R1_001.fastq

seqtk sample -s999 iCELL8_01_S1_L001_R1_001.fastq.paired.fq 250000 > iCELL8_01_S1_L001_R1_001.fastq         
seqtk sample -s100 iCELL8_01_S1_L001_R2_001_2_chr21_R2.fastq.paired.fq 250000 > iCELL8_01_S1_L002_R2_001.fastq
seqtk sample -s100 iCELL8_01_S1_L001_R1_001.fastq.paired.fq 250000 > iCELL8_01_S1_L002_R1_001.fastq
seqtk sample -s999 iCELL8_01_S1_L001_R2_001_2_chr21_R2.fastq.paired.fq 250000 > iCELL8_01_S1_L001_R2_001.fastq
