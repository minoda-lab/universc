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
#Script "MakeNewReference.pl" given a genome file (FASTA), annotation file (GTF), and a genome name, generates a reference for cellranger.
###########



#####Options#####
print "\n\#\#\#\#\#Running MakeNewReference\#\#\#\#\#\n";
print "checking options.\n";


#setting the default values.
my $genome = "";
my $annotation = "";
my $name = "";

#making the options into external arguments.
GetOptions (
        'genome=s' => \$genome,
        'annotation=s' => \$annotation,
        'name=s' => \$name
        );

#checking for sequence type
if (!$genome) {
        die "USAGE: option --genome <FASTA file> is required.\n";
}
elsif (!$annotation){
        die "USAGE: option --annotation <GTF file> is required.\n";
}
elsif (!$name) {
        die "USAGE: option --name <reference name> is required.\n";
}

print "check complete.\n\n";
##########



#####MAIN#####
print "making cellranger reference.\n";

system "cellranger mkref --genome=$name --fasta=$genome --genes=$annotation";

print "reference made.\n";
print "\#\#\#\#\#MakeNewReference COMPLETE\#\#\#\#\#\n\n";
##########



#####SUB#####
##########
__END__
