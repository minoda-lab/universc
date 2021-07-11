#!/bin/bash

PWD=$(shell pwd)

all: test-data

test-data:
	@echo "generating test data from public repository"
	@if [ ! -f Smartseq3_Fibroblasts_GelCut_small_R1.fastq.gz ]; then\
		bash prepare_files.sh
	else\
		echo "test data pre-generated";\
	fi
