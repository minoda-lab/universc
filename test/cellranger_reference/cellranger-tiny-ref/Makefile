#!/bin/bash

PWD=$(shell pwd)

all: reference star

reference:
	@echo "generating cellranger references"
	@if [ ! -d 1.2.0 ]; then\
		echo cellranger mkref --genome=1.2.0 --fasta=genome-1.2.0.fa --genes=genes-1.2.0.gtf;\
		rm -rf genome-1.2.0.fa genes-1.2.0.gtf;\
	else\
		echo "reference 1.2.0 pre-generated";\
	fi
	@if [ ! -d 3.0.0 ]; then\
		echo cellranger mkref --genome=3.0.0 --fasta=genome-3.0.0.fa --genes=genes-3.0.0.gtf;\
		rm -rf genome-3.0.0.fa genes-3.0.0.gtf;\
	else\
		echo "reference 3.0.0 pre-generated";\
	fi

star:
        @echo "check for SA files"
        # used to export to PATH for testing on SGE server
        @cellrangerversion=`cellranger count --version | head -n 2 | tail -n 1 | cut -f2 -d'(' | cut -f1 -d')'`
        @cellrangerpath=`which cellranger`
        # set up cellranger reference
        @if [[ ! -f ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/3.0.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ]]; then
                rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/3.0.0/star/SA ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/3.0.0/star/SA
        fi
        @if [[ ! -f ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/1.2.0/star/SA ]] && [[ -f $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ]]; then
                rsync $(dirname $cellrangerpath)/cellranger-tiny-ref/1.2.0/star/SA ${cellrangerpath}-${cellrangerversion}/cellranger-tiny-ref/1.2.0/star/SA
        fi

##########################################
# Generic rules
##########################################

clean:
	echo "moving fasta and gtf from cellranger reference"
	mv 1.2.0/fasta/genome.fa genome-1.2.0.fa
	mv 1.2.0/genes/genes.gtf genes-1.2.0.gtf
	@if [ -d 3.0.0 ]; then\
		rm -rf 1.2.0;\
	fi
	mv 3.0.0/fasta/genome.fa genome-3.0.0.fa
	mv 3.0.0/genes/genes.gtf genes-3.0.0.gtf
	@if [ -d 3.0.0 ]; then\
		rm -rf 3.0.0;\
	fi	

mrproper: clean
	rm -f *~ $(EXEC)
