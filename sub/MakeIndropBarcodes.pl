#!/usr/bin/perl

#####################################
#                                   #
#     Written by Kai Battenberg     #
#     Plant Sciences UC Davis       #
#                                   #
#####################################

use strict;
use warnings;

#####Options#####
my $gel_list1 = $ARGV[0];
my $gel_list2 = $ARGV[1];
my $version = $ARGV[2];
my $out = $ARGV[3];
##########



#####MAIN#####
my $list1 = `cat $gel_list1`;
chomp $list1;
my @list1 = split (/\n/, $list1);

my $list2 = `cat $gel_list2`;
chomp $list2;
my @list2 = split (/\n/, $list2);

my $length = scalar (@list1);



my @list12;
if ($version =~ m/v2/) {
    for (my $i = 0; $i < $length; $i++) {
        my $element1 = substr ($list1[$i], 0, 8);
        my $element12 = $element1.$list2[$i];
        push (@list12, $element12);
    }
}
elsif ($version =~ m/v3/) {
    foreach my $element1 (@list1) {
        $element1 = substr ($element1, 0, 8);
        foreach my $element2 (@list2) {
            my $element12 = $element1.$element2;
            push (@list12, $element12);
        }
    }
}

my %list12;
foreach my $element12 (sort @list12) {
    if (!exists $list12{$element12}) {
        $list12{$element12} = 1;
    }
    else {
        $list12{$element12}++
    }
}

my $outfile = $out."/inDrop-".$version."_barcodes.txt";
open (OUT, ">$outfile") or die "cannot open $outfile.\n";
foreach my $element12 (sort keys %list12) {
    if ($list12{$element12} < 2) {
        print OUT "$element12\n";
    }
}
close (OUT);
##########
__END__
