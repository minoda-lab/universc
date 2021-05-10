#!/usr/bin/perl

#########################################
#                                       #
#     Written by Kai Battenberg         #
#     Plant Symbiosis Research Team     #
#                                       #
#########################################

use strict;
use warnings;

#####Options#####
my $barcodefile = $ARGV[0];
my $adjust_length = $ARGV[1];
my $cellranger_dir = $ARGV[2];
##########



#####Main#####
#make a hash of original and adjusted barcodes
my %barcodes;
open (ORIGINAL, "<$barcodefile") or die "cannot open $barcodefile.\n";
while (my $line = <ORIGINAL>) {
    my $original = $line;
    $original =~ s/\r//sig;
    $original =~ s/\n//sig;
    
    my $adjusted = $original;
    if ($adjust_length >= 0) {
        my $length = $adjust_length;
        $adjusted = substr ($adjusted, $length);
    }
    else {
        my $length = abs ($adjust_length);
        my $As = "A" x $length;
        $adjusted = $As.$adjusted;
    }
    
    $barcodes{$adjusted} = $original;
}
close (ORIGINAL);

#get modified barcode list (raw)
my $raw_file = $cellranger_dir."/outs/raw_feature_bc_matrix/barcodes.tsv.gz";
my $modified_barcodes_raw = `zcat $raw_file | cut -f1 -d'-'`;
$modified_barcodes_raw =~ s/\n$//;
my @modified_barcodes_raw = split (/\n/, $modified_barcodes_raw);

#get modified barcode list (filtered)
my $filtered_file = $cellranger_dir."/outs/filtered_feature_bc_matrix/barcodes.tsv.gz";
my $modified_barcodes_filtered = `zcat $filtered_file | cut -f1 -d'-'`;
$modified_barcodes_filtered =~ s/\n$//;
my @modified_barcodes_filtered = split (/\n/, $modified_barcodes_filtered);

#replace barcodes
my $temp = "barcodes.txt";
open (TEMP, ">$temp") or die "cannot open $temp.\n";
foreach my $barcode (@modified_barcodes_raw) {
    print TEMP "$barcodes{$barcode}-1\n";
}
close (TEMP);
system "gzip -c $temp > $raw_file";

open (TEMP, ">$temp") or die "cannot open $temp.\n";
foreach my $barcode (@modified_barcodes_filtered) {
    print TEMP "$barcodes{$barcode}-1\n";
}
close (TEMP);
system "gzip -c $temp > $filtered_file";
unlink ($temp);
##########
