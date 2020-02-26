#!/bin/bash

PWD=$(shell pwd)

all: $(EXEC)

louvain:
	echo "generating cellranger references"
	cellranger mkref --genome=1.2.0 --fasta=genome-1.2.0.fa --genes=genes-1.2.0.gtf
	rm -rf genome-1.2.0.fa genes-1.2.0.gtf
	cellranger mkref --genome=3.0.0 --fasta=genome-3.0.0.fa --genes=genes-3.0.0.gtf
	rm -rf genome-3.0.0.fa genes-3.0.0.gtf

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	echo "moving fasta and gtf from cellranger reference"
	mv 1.2.0/fasta/genome.fa genome-1.2.0.fa
	mv 1.2.0/genes/genes.gtf genes-1.2.0.gtf
	mv 3.0.0/fasta/genome.fa genome-3.0.0.fa
	mv 3.0.0/genes/genes.gtf genes-3.0.0.gtf

mrproper: clean
	rm -f *~ $(EXEC)
