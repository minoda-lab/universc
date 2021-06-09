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
#Script "FilterSmartSeqReadUMI.pl" given a set of I1, I2, R1, R2, tag, will generate
#   1. R1 file parsed and converted based on tag and TSO sequence respectively
#   2. R2, I1, and I2 files with reads corresponding to R1
###########



#####Options#####
#target tag
my $tag = "";

#target tso
my $tso_10x = "TTTCTTATATGGG";
my $tso_10x_q = "IIIIIIIIIIIII";

#setting the default values.
my $i1 = "";
my $i2 = "";
my $r1 = "";
my $r2 = "";
my $i1_out = "parsed_I1.fastq";
my $i2_out = "parsed_I2.fastq";
my $r1_out = "parsed_R1.fastq";
my $r2_out = "parsed_R2.fastq";
my $out_dir = "";

#making the options into external arguments.
GetOptions (
	'tag=s' => \$tag,
	'i1=s' => \$i1,
	'i2=s' => \$i2,
	'r1=s' => \$r1,
	'r2=s' => \$r2,
	'out_dir=s' => \$out_dir
	);

#Checking options
if (!$i1) {
	die "USAGE: option --i1 <I1.fastq> is required.\n";
}
elsif (!$i2) {
	die "USAGE: option --i2 <I2.fastq> is required.\n";
}
elsif (!$r1) {
	die "USAGE: option --r1 <R1.fastq> is required.\n";
}
elsif (!$r2) {
	die "USAGE: option --r2 <R2.fastq> is required.\n";
}
elsif (!$tag) {
	die "USAGE: option --tag <TARGET_TAG_SEQUENCE> is required.\n";
}
elsif (!$out_dir) {
	die "USAGE: option --out_dir <OUTPUT_DIRECTORY> is required.\n";
}

#set output file names
$i1_out = $out_dir."/".$i1_out;
$i1_out =~ s/\/\//\//;
$i2_out = $out_dir."/".$i2_out;
$i2_out =~ s/\/\//\//;
$r1_out = $out_dir."/".$r1_out;
$r1_out =~ s/\/\//\//;
$r2_out = $out_dir."/".$r2_out;
$r2_out =~ s/\/\//\//;

#set technology
my $technology = "";
if ($tag eq "AAGCAGTGGTATCAACGCAGAGTAC") {
	$technology = "SmartSeq2";
}
elsif ($tag eq "ATTGCGCAATG") {
	$technology = "SmartSeq3";
}
##########



#####MAIN#####
#parse reads to keep
my $readcount = 0;
open (I1, "<", $i1) or die "cannot open $i1.\n";
open (I2, "<", $i2) or die "cannot open $i2.\n";
open (R1, "<", $r1) or die "cannot open $r1.\n";
open (R2, "<", $r2) or die "cannot open $r2.\n";
open (I1OUT, ">", $i1_out) or die "cannot open $i1_out.\n";
open (I2OUT, ">", $i2_out) or die "cannot open $i2_out.\n";
open (R1OUT, ">", $r1_out) or die "cannot open $r1_out.\n";
open (R2OUT, ">", $r2_out) or die "cannot open $r2_out.\n";
while (my $line = <R1>) {
	$readcount++;
	
	#read in data
	my $i1_head = <I1>;
	my $i1_seq = <I1>;
	my $i1_plus = <I1>;
	my $i1_q = <I1>;
	my $i2_head = <I2>;
	my $i2_seq = <I2>;
	my $i2_plus = <I2>;
	my $i2_q = <I2>;
	my $r1_head = $line;
	my $r1_seq = <R1>;
	my $r1_plus = <R1>;
	my $r1_q = <R1>;
	my $r2_head = <R2>;
	my $r2_seq = <R2>;
	my $r2_plus = <R2>;
	my $r2_q = <R2>;
	
	#remove reads without a tag
	if ($r1_seq !~ m/$tag/) {
		next;
	}
	
	#trim r1 data by tag
	my $r1_trim_seq = $r1_seq;
	my @trim = split (/$tag/, $r1_trim_seq);
	shift @trim;
	$r1_trim_seq = join ("$tag", @trim);
	$r1_trim_seq = $r1_trim_seq;
	my $r1_trim_q = $r1_q;
	chomp $r1_trim_q;
	$r1_trim_q = reverse $r1_trim_q;
	$r1_trim_q = substr ($r1_trim_q, 0, length($r1_trim_seq) - 1);
	$r1_trim_q = reverse $r1_trim_q;
	$r1_trim_q = $r1_trim_q."\n";
	
	#remove reads too short after trimming
	if ($technology eq "SmartSeq2" && length($r1_trim_seq) < 3) {
		next;
	}
	elsif ($technology eq "SmartSeq3" && length($r1_trim_seq) < 11) {
		next;
	}
	
	#replacing tso
	my $r1_trim_swop_seq;
	my $r1_trim_swop_q;
	if ($technology eq "SmartSeq2") {
		$r1_trim_swop_seq = $tso_10x.(substr ($r1_trim_seq, 3));
		$r1_trim_swop_q = $tso_10x_q.(substr ($r1_trim_q, 3));
	}
	elsif ($technology eq "SmartSeq3") {
		$r1_trim_swop_seq = (substr ($r1_trim_seq, 0, 8)).$tso_10x.(substr ($r1_trim_seq, 11));
		$r1_trim_swop_q = (substr ($r1_trim_q, 0, 8)).$tso_10x_q.(substr ($r1_trim_q, 11));
	}
	
	#print out data
	print I1OUT "$i1_head";
	print I1OUT "$i1_seq";
	print I1OUT "$i1_plus";
	print I1OUT "$i1_q";
	print I2OUT "$i2_head";
	print I2OUT "$i2_seq";
	print I2OUT "$i2_plus";
	print I2OUT "$i2_q";
	print R1OUT "$r1_head";
	print R1OUT "$r1_trim_swop_seq";
	print R1OUT "$r1_plus";
	print R1OUT "$r1_trim_swop_q";
	print R2OUT "$r2_head";
	print R2OUT "$r2_seq";
	print R2OUT "$r2_plus";
	print R2OUT "$r2_q";
}
close (I1);
close (I2);
close (R1);
close (R2);
close (I1OUT);
close (I2OUT);
close (R1OUT);
close (R2OUT);
##########
