wget https://www.ncbi.nlm.nih.gov/geo/download/\?acc\=GSM1629192\&format\=file\&file\=GSM1629192%5FPure%5FHumanMouse%2Ebam
mv index.html\?acc=GSM1629192\&format=file\&file=GSM1629192%5FPure%5FHumanMouse%2Ebam GSM162919.bam
samtools sort -n GSM162919.bam > GSM162919.qsort
samtools view  GSM162919.qsort  HUMAN_21:9825832-48085036 > GSM162919.qsort2
samtools sort -O BAM GSM162919.bam > GSM162919.sort.bam
samtools index GSM162919.sort.bam
samtools view  GSM162919.sort.bam  HUMAN_21:9825832-48085036 > GSM162919.chr21.bam
samtools view -O BAM  GSM162919.sort.bam  HUMAN_21:9825832-48085036 > GSM162919.chr21.sort.bam
samtools sort -n GSM162919.chr21.sort.bam -o GSM162919.chr21.qsort.bam
bedtools bamtofastq -i GSM162919.chr21.qsort.bam -fq GSM1629192_chr21_R1.fastq
mv GSM1629192_chr21_R1.fastq GSM1629192_chr21_R2.fastq
fastq-dump -F --split-files SRR1873277
fastq_pair GSM1629192_chr21_R2.fastq SRR1873277_1.fastq
head -n 117060 SRR1873277_1.fastq.paired.fq 117060 > SRR1873277_1.fastq.paired.fq
head -n 117060 GSM1629192_chr21_R2.fastq.paired.fq > GSM1629192_chr21_R2.fastq.paired.fq
cp SRR1873277_1.fastq.paired.fq  GSM1629192_chr21_R2.fastq.paired.fq ~/repos/universc/test/shared/dropseq-test
cp SRR1873277_1.fastq.paired.fq  GSM1629192_chr21_R2.fastq.paired.fq ~/repos/universc/test/shared/dropseq-test
mv SRR1873277_1.fastq.paired.fq SRR1873277_R1.fastq
mv GSM1629192_chr21_R2.fastq.paired.fq  SRR1873277_R2.fastq
mv GSM1629192_chr21_R2.fastq.paired.fq  SRR1873277_R2.fastq
