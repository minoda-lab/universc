#!/usr/bin/perl

#########################################
#                                       #
#     Written by Kai Battenberg         #
#     Plant Symbiosis Research Team     #
#                                       #
#########################################

use strict;
use warnings;
use List::Util qw( min max sum);

#####Options#####
my $barcodefile=$ARGV[0];
my $adjust_length=$ARGV[1];
my $indir = $ARGV[2];
##########



#####Check input#####
#collect files needed
my $feature_file = $indir."/outs/raw_feature_bc_matrix/features.tsv.gz";
my $raw_barcode_file = $indir."/outs/raw_feature_bc_matrix/barcodes.tsv.gz";
my $filtered_barcode_file = $indir."/outs/filtered_feature_bc_matrix/barcodes.tsv.gz";
my $matrix_file = $indir."/outs/filtered_feature_bc_matrix/matrix.mtx.gz";
my $summary_file = $indir."/outs/metrics_summary.csv";
my $bam_file = $indir."/outs/possorted_genome_bam.bam";
my $out_file = $indir."/outs/basic_stats.txt";

if (! -e $indir) {
    die "Error: Data directory needs to be specified.\n";
}
elsif (! -e $feature_file || ! -e $raw_barcode_file || ! -e $filtered_barcode_file || ! -e $matrix_file || ! -e $summary_file || ! -e $bam_file) {
    die "Error: Iether selected directory is not a cellranger output or some critical files are missing.\n";
}
##########



#####MAIN#####
#make a hash of original and adjusted barcodes
my %adjusted2original;
open (ORIGINAL, "<$barcodefile") or die "cannot open $barcodefile.\n";
while (my $line = <ORIGINAL>) {
    my $original = $line;
    $original =~ s/\r//sig;
    $original =~ s/\n//sig;
    
    my $adjusted = $original;
    if ($adjust_length >= 0) {
        my $length = $adjust_length;
        $adjusted = substr ($adjusted, 0, $length);
    }
    else {
        my $length = abs ($adjust_length);
        my $As = "A" x $length;
        $adjusted = $As.$adjusted;
    }
    
    $adjusted2original{$adjusted} = $original;
}
close (ORIGINAL);

#gene count
print " getting gene count.\n";

my $ref_gene_count = `zcat $feature_file | wc -l`;
chomp $ref_gene_count;

#cell count
print " getting cell count.\n";

my $cellcount_raw;
$cellcount_raw = `zcat $raw_barcode_file | wc -l`;
chomp $cellcount_raw;

my $cellcount_filtered;
$cellcount_filtered = `zcat $filtered_barcode_file | wc -l`;
chomp $cellcount_filtered;

my %percell;
my $barcodes;
$barcodes = `zcat $filtered_barcode_file | cut -f1 -d'-'`;
chomp $barcodes;
my @barcodes = split (/\n/, $barcodes);

foreach my $barcode (@barcodes) {
    my %barcode;
    $percell{$barcode} = \%barcode;
}

#raw read count
print " getting raw read count.\n";

my $readcount_raw = `tail -n 1 $summary_file`;
chomp $readcount_raw;

my @readcount_raw = split (//, $readcount_raw);
my $quotestat = 0;
foreach my $character (@readcount_raw) {
	if ($character =~ m/\"/ && $quotestat == 0) {
		$character = "";
		$quotestat = 1;
		next;
	}
	elsif ($character =~ m/\"/ && $quotestat == 1) {
		$character = "";
		$quotestat = 0;
		next;
	}
	elsif ($character =~ m/\,/ && $quotestat == 1) {
		$character = "";
		next;
	}
	else {
		next;
	}
}
$readcount_raw = join ("", @readcount_raw);
$readcount_raw = (split (/\,/, $readcount_raw))[3];

#assigned read count
print " getting assigned read count.\n";

my $assigned_reads_percell = `samtools view $bam_file | grep 'CB:Z:' | sed 's\/CB:Z:\\([ACGT]*\\).*/\\1/' | awk '{print \$1, \$NF}' | sort | uniq | cut -f2 -d' ' | sort | uniq -c`;
chomp $assigned_reads_percell;
my @assigned_reads_percell = split (/\n/, $assigned_reads_percell);

foreach my $barcode (@assigned_reads_percell) {
    my $count = $barcode;
    $count =~ s/ +/ /g;
    $count =~ s/^ //;
    $count =~ s/ $//;
    my @count = split (/ /, $count);
    
    if (exists $percell{ $adjusted2original{$count[1]} }) {
        my %count;
        $count{ASSIGNED} = $count[0];
        $percell{ $adjusted2original{$count[1]} } = \%count;
    }
}
foreach my $barcode (sort keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    if (!$barcode{ASSIGNED}) {
        $barcode{ASSIGNED} = 0;
    }
    $percell{$barcode} = \%barcode;
}

#mapped read count
print " getting mapped read count.\n";

my $mapped_reads_percell = `samtools view -F 4 $bam_file | grep 'CB:Z:' | sed 's\/CB:Z:\\([ACGT]*\\).*/\\1/' | awk '{print \$1, \$NF}' | sort | uniq | cut -f2 -d' ' | sort | uniq -c`;
chomp $mapped_reads_percell;
my @mapped_reads_percell = split (/\n/, $mapped_reads_percell);

foreach my $barcode (@mapped_reads_percell) {
    my $count = $barcode;
    $count =~ s/ +/ /g;
    $count =~ s/^ //;
    $count =~ s/ $//;
    my @count = split (/ /, $count);

    if (exists $percell{ $adjusted2original{$count[1]} }) {
        my %count = %{ $percell{ $adjusted2original{$count[1]} } };
        $count{MAPPED} = $count[0];
        $percell{ $adjusted2original{$count[1]} } = \%count;
    }
}
foreach my $barcode (sort keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    if (!$barcode{MAPPED}) {
        $barcode{MAPPED} = 0;
    }
    $percell{$barcode} = \%barcode;
}

#UMI count
print " getting mapped UMI count.\n";

my %barcodes;
my $id = 0;
foreach my $barcode (@barcodes) {
    $id++;
    $barcodes{$id} = $barcode;
}

my $umicount = `zcat $matrix_file | tail -n +4 | cut -f2,3 -d ' '`;
chomp $umicount;
my @umicount = split (/\n/, $umicount);
my %umicount;
foreach my $entry (@umicount) {
    my @data = split (/ /, $entry);
    if (!exists $umicount{$data[0]}) {
        $umicount{$data[0]} = $data[1];
    }
    else {
        $umicount{$data[0]} = $umicount{$data[0]} + $data[1];
    }
}
foreach my $id (keys %umicount) {
    if (exists $percell{ $barcodes{$id} }) {
        my %data = %{ $percell{ $barcodes{$id} } };
        $data{UMI} = $umicount{$id};
        
        $percell{ $barcodes{$id} } = \%data;
    }
}
foreach my $barcode (sort keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    if (!$barcode{UMI}) {
        $barcode{UMI} = 0;
    }
    $percell{$barcode} = \%barcode;
}

#Gene count per cell
print " getting gene count per cell.\n";

my @genecount = @umicount;
my %genecount;
foreach my $entry (@genecount) {
    my @data = split (/ /, $entry);
    if (!exists $genecount{$data[0]}) {
        $genecount{$data[0]} = 1;
    }
    else {
        $genecount{$data[0]}++;
    }
}
foreach my $id (keys %genecount) {
    if (exists $percell{ $barcodes{$id} }) {
        my %data = %{ $percell{ $barcodes{$id} } };
        $data{GENE} = $genecount{$id};
        
        $percell{ $barcodes{$id} } = \%data;
    }
}
foreach my $barcode (sort keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    if (!$barcode{GENE}) {
        $barcode{GENE} = 0;
    }
    $percell{$barcode} = \%barcode;
}

#Total gene count
print " getting total gene count.\n";

my $genes_total = `zcat $matrix_file | tail -n +4 | cut -f1 -d ' ' | sort | uniq | wc -l`;
chomp $genes_total;

#Total, means, and medinans
print " getting total, mean, and medians.\n";

my @reads_assigned;
my @reads_mapped;
my @umis;
my @genes;
foreach my $barcode (keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    push (@reads_assigned, $barcode{ASSIGNED});
    push (@reads_mapped, $barcode{MAPPED});
    push (@umis, $barcode{UMI});
    push (@genes, $barcode{GENE});
}

my $reads_assigned_total = sum (@reads_assigned);
my $reads_assigned_mean = &mean (@reads_assigned);
my $reads_assigned_median = &median (@reads_assigned);

my $reads_mapped_total = sum (@reads_mapped);
my $reads_mapped_mean = &mean (@reads_mapped);
my $reads_mapped_median = &median (@reads_mapped);

my $umis_total = sum (@umis);
my $umis_mean = &mean (@umis);
my $umis_median = &median (@umis);

my $genes_mean = &mean (@genes);
my $genes_median = &median (@genes);

#getting percentages
print " getting percentages.\n";

my $percent_assigned = sprintf("%.4f", $reads_assigned_total/$readcount_raw)*100;
$percent_assigned = $percent_assigned."\%";
my $percent_mapped = sprintf("%.4f", $reads_mapped_total/$reads_assigned_total)*100;
$percent_mapped = $percent_mapped."\%";
my $reads_per_umi = sprintf("%.2f", $reads_mapped_total/$umis_total);

#output
print " printing out summary.\n";

open (OUT, ">$out_file") or die "cannot open $out_file.\n";
print OUT "\#\#\#\#\#BASIC STATS\#\#\#\#\#\n";
print OUT "Reference gene count:\t$ref_gene_count\n";
print OUT "Cell count:\t(raw)$cellcount_raw\t(filtered)$cellcount_filtered\n";
print OUT "Read count(raw):\t$readcount_raw\n";
print OUT "Read count(assigned):\t(total)$reads_assigned_total\t(per cell mean)$reads_assigned_mean\t(per cell median)$reads_assigned_median\n";
print OUT "Read count(mapped):\t(total)$reads_mapped_total\t(per cell mean)$reads_mapped_mean\t(per cell median)$reads_mapped_median\n";
print OUT "UMI count(mapped to gene):\t(total)$umis_total\t(per cell mean)$umis_mean\t(per cell median)$umis_median\n";
print OUT "Gene count:\t(total)$genes_total\t(per cell mean)$genes_mean\t(per cell median)$genes_median\n";
print OUT "\%Reads assigned/raw:\t$percent_assigned\n";
print OUT "\%Reads mapped/assigned:\t$percent_mapped\n";
print OUT "reads/UMI:\t$reads_per_umi\n";
print OUT "\n";
print OUT "\#\#\#\#\#PER CELL\#\#\#\#\#\n";
print OUT "barcode\tRead count (assigned)\tRead count (mapped)\tUMI count(mapped)\tGene count\n";
foreach my $barcode (sort keys %percell) {
    my %barcode = %{ $percell{$barcode} };
    print OUT "$barcode\t$barcode{ASSIGNED}\t$barcode{MAPPED}\t$barcode{UMI}\t$barcode{GENE}\n";
}
close (OUT);
##########



#####SUB#####
sub mean {
    my $raw = sum(@_)/@_;
    my $round = sprintf("%.1f", $raw); 
    return $round;
}

sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) {
        return $vals[int($len/2)];
    }
    else {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
##########
__END__
