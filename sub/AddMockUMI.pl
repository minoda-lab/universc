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
#Script "add_mock_UMI.pl" given a FASTQ file, starting position, and desired UMI length, generates a new FASTQ file with a mock UMI inserted.
###########



######Options#####
my $fastq_in = "";
my $fastq_out = "mock_UMI.fastq";
my $head_length = ""; #number of characters to have before insertion of mock UMI.
my $umi_length = ""; #length of desired UMI
##########



#####Checking Options#####
#making the options into external arguments.
GetOptions (
	'fastq=s' => \$fastq_in,
	'head_length=s' => \$head_length,
	'umi_length=s' => \$umi_length
	);

#checking for required options.
if (!$fastq_in) {
	die "USAGE: option --fastq <FASTQ FILE> is required.\n";
}
elsif (!$head_length) {
	die "USAGE: option --head_length <INTEGER> is required.\n";
}
elsif (!$umi_length) {
	die "USAGE: option --umi_length <INTEGER> is required.\n";
}

#checking head length
if ( $head_length !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --head_length needs to be numeric\n";
}
elsif ( $head_length < 0 || $head_length !~ /^-?\d+$/ ) {
	die "Error: option --head_length needs to be an integer 0 or greater\n";
}

#checking umi length
if ( $umi_length !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --umi_length needs to be numeric\n";
}
elsif ( $head_length < 1 || $head_length !~ /^-?\d+$/ ) {
	die "Error: option --umi_length needs to be an integer 0 or greater\n";
}
##########



#####Main#####
my %nucleotide = (
	'00' => "A",
	'01' => "C",
	'10' => "G",
	'11' => "T"
);

my $insert_q = "I" x $umi_length;
my $read_count = -1;
open (IN, "<", $fastq_in) or die "cannot open $fastq_in.\n";
open (OUT, ">", $fastq_out) or die "cannot open $fastq_out.\n";
while (my $line = <IN>) {
	$read_count++;
	
	#get binary read count
	my $binary = sprintf ("%b", $read_count);
	$binary = (0 x (2 * $umi_length)).$binary;
	$binary = reverse(substr (reverse ($binary), 0, 2 * $umi_length));
	my @binary = ( $binary =~ m/../g );
	
	#get mock UMI
	my $mock_umi = "";
	foreach my $character (@binary) {
		$mock_umi = $mock_umi.$nucleotide{$character};
	}
	
	#get data from FASTQ file
	my $header = $line;
	chomp $header;
	my $seq = <IN>;
	chomp $seq;
	my $plus = <IN>;
	my $q = <IN>;
	chomp $q;
	
	#insert 
	my $seq_head = substr($seq, 0, $head_length);
	my $seq_tail = reverse(substr (reverse ($seq), 0, length($seq) - $head_length));
	my $q_head = substr($q, 0, $head_length);
	my $q_tail = reverse(substr (reverse ($q), 0, length($q) - $head_length));
	
	print OUT "$header\n";
	print OUT "$seq_head$mock_umi$seq_tail\n";
	print OUT "$plus";
	print OUT "$q_head$insert_q$q_tail\n";
}
close (IN);
close (OUT);
##########
__END__
