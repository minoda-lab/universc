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
#Script "FilterSmartSeqReadUMI.pl" given a set of I1, I2, R1, and R2, trims and parses R1 and only keep the corresponding reads in the other files.
###########



#####Options#####
#SmartSeq3 tag
my $ss3t = "ATTGCGCAATG";

#setting the default values.
my $i1 = "";
my $i2 = "";
my $r1 = "";
my $r2 = "";
my $i1_out = "SmartSeq3_parsed_I1.fastq";
my $i2_out = "SmartSeq3_parsed_I2.fastq";
my $r1_out = "SmartSeq3_parsed_R1.fastq";
my $r2_out = "SmartSeq3_parsed_R2.fastq";
my $out_dir = "";

#making the options into external arguments.
GetOptions (
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
##########



#####MAIN#####
#get list of reads to keep
my @keepreads = map { ($_ + 2) / 4 } split(/\n/, `grep $ss3t $r1 -n | cut -f1 -d\:`);
my %keepreads = map { $_ => 1 } @keepreads;

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
	
	if (!exists $keepreads{$readcount}) {
		next;
	}
	else {
		#trim r1 data
		my @trim = split (/$ss3t/, $r1_seq);
		shift @trim;
		my $new_r1_seq = join ("$ss3t", @trim);
		my $new_r1_q = $r1_q;
		chomp $new_r1_q;
		$new_r1_q = reverse $new_r1_q;
		$new_r1_q = substr ($new_r1_q, 0, length($new_r1_seq) - 1);
		$new_r1_q = reverse $new_r1_q;
		$new_r1_q = $new_r1_q."\n";
		
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
		print R1OUT "$new_r1_seq";
		print R1OUT "$r1_plus";
		print R1OUT "$new_r1_q";
		print R2OUT "$r2_head";
		print R2OUT "$r2_seq";
		print R2OUT "$r2_plus";
		print R2OUT "$r2_q";
	}
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
