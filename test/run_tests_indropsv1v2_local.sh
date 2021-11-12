#! /bin/bash

bash launch_universc.sh --id "test-in-v2" --technology "indrop2" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "/home/tom/datasets/indropsv2_public_data/lib1_L001" \
 "/home/tom/datasets/indropsv2_public_data/lib1_L002" "/home/tom/datasets/indropsv2_public_data/lib1_L003" "/home/tom/datasets/indropsv2_public_data/lib1_L004" \
 --per-cell-data --jobmode "local" --localcores 8

bash launch_universc.sh --id "test-in-v1" --technology "indrop1" \
 --reference "test/cellranger_reference/cellranger-tiny-ref/3.0.0" \
 --file "/home/tom/datasets/indropsv1_public_data/sample_S1_L001" \
 --per-cell-data --jobmode "local" --localcores 8

rename  "s/S1_L001/L001/" /home/tom/datasets/indropsv2_public_data/sample*fastq
rename  "s/S2_L002/L002/" /home/tom/datasets/indropsv2_public_data/sample*fastq
rename  "s/S3_L003/L003/" /home/tom/datasets/indropsv2_public_data/sample*fastq
rename  "s/S4_L004/L004/" /home/tom/datasets/indropsv2_public_data/sample*fastq
