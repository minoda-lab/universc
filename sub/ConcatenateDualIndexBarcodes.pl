#!/usr/bin/perl

#########################################
#                                       #
#     Written by Kai Battenberg         #
#     Plant Symbiosis Research Team     #
#                                       #
#########################################

use strict;
use warnings;
use Getopt::Long;

#####SCRIPT DESCRIPTION#####
#Script "ConcatenateDualIndexBarcodes.pl" given a reference file, and a list of fastq files, concatenates all the additive sequences to add in front of each reference read.
###########



#####Options#####
#setting the default values.
my $ref_fastq = "";
my @additives = ();
my $concatenated_fastq = "Concatenated_File.fastq";
my $out_dir = "";

#making the options into external arguments.
GetOptions (
	'ref_fastq=s' => \$ref_fastq,
	'additive=s' => \@additives,
	'out_dir=s' => \$out_dir
	);

#Checking options
if (!$ref_fastq) {
	die "USAGE: option --ref_fastq <Reference.fastq> is required.\n";
}
elsif (scalar @additives < 1) {
	die "USAGE: option --additive <additional.fastq> is required.\n";
}
elsif (!$out_dir) {
	die "USAGE: option --out_dir <OUTPUT_DIRECTORY> is required.\n";
}

$concatenated_fastq = $out_dir."/".$concatenated_fastq;
$concatenated_fastq =~ s/\/\//\//;
##########



#####MAIN#####
#read in additive data
my %additives_seq;
my %additives_q;
foreach my $file (@additives) {
	my $readcount = 0;
	
	open (ADD, "<", $file) or die "cannot open $file.\n";
	while (my $line = <ADD>) {
		$readcount++;
		
		my $head = $line;
		my $seq = <ADD>;
		my $plus = <ADD>;
		my $q = <ADD>;
		
		chomp $seq;
		chomp $q;
		
		#save additive data
		if (!$additives_seq{$readcount}) {
			$additives_seq{$readcount} = $seq;
		}
		else {
			$additives_seq{$readcount} = $additives_seq{$readcount}.$seq;
		}
		
		if (!$additives_q{$readcount}) {
			$additives_q{$readcount} = $q;
		}
		else {
			$additives_q{$readcount} = $additives_q{$readcount}.$q;
		}
	}
	close (ADD);
}

#add data to reference
my $ref_readcount = 0;
open (REF, "<", $ref_fastq) or die "cannot open $ref_fastq.\n";
open (CAT, ">", $concatenated_fastq) or die "cannot open $concatenated_fastq.\n";
while (my $line = <REF>) {
	$ref_readcount++;
	
	my $head = $line;
	my $seq = <REF>;
	my $plus = <REF>;
	my $q = <REF>;
	
	my $new_seq = $additives_seq{$ref_readcount}.$seq;
	my $new_q = $additives_q{$ref_readcount}.$q;
	
	print CAT "$head";
	print CAT "$new_seq";
	print CAT "$plus";
	print CAT "$new_q";
}
close (REF);
close (CAT);
##########
