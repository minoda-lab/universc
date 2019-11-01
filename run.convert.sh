#!bin/bash

export PATH=/home/kai.b/tools/cellranger-3.1.0:$PATH

bash /home/kai.b/tools/launch_universc.sh \
	-R1 /home/kai.b/scRNAseq_ICELL8_3primeNoUMI/TEST/01_INPUT/I03Ljmini_S1_L002_R1_001.fastq \
	-R2 /home/kai.b/scRNAseq_ICELL8_3primeNoUMI/TEST/01_INPUT/I03Ljmini_S1_L002_R2_001.fastq \
	-s \
	-t icell8 \
	-i I03Ljmini \
	-r /home/kai.b/tools/cellranger_references/Ljaponicus_v2.5 \
	-c SC3Pv2 \
	-j sge \
	-w convert \
	--verbose
