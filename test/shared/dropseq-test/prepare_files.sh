wget https://www.ncbi.nlm.nih.gov/geo/download/\?acc\=GSM1629192\&format\=file\&file\=GSM1629192%5FPure%5FHumanMouse%2Ebam
mv index.html\?acc=GSM1629192\&format=file\&file=GSM1629192%5FPure%5FHumanMouse%2Ebam GSM1629192.bam
samtools sort -n GSM1629192.bam > GSM1629192.qsort
samtools view  GSM1629192.qsort  HUMAN_21:9825832-48085036 > GSM1629192.qsort2
samtools sort -O BAM GSM1629192.bam > GSM1629192.sort.bam
samtools index GSM1629192.sort.bam
samtools view  GSM1629192.sort.bam  HUMAN_21:9825832-48085036 > GSM1629192.chr21.bam
samtools view -O BAM  GSM1629192.sort.bam  HUMAN_21:9825832-48085036 > GSM1629192.chr21.sort.bam
samtools sort -n GSM1629192.chr21.sort.bam -o GSM1629192.chr21.qsort.bam
bedtools bamtofastq -i GSM1629192.chr21.qsort.bam -fq GSM1629192_chr21_R1.fastq
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
