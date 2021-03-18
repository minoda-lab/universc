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
#Script "TrimSeq4scRNAseq.pl" trims the R2 file and pairs it with the untrimmed R1 file.
###########



######Options#####
print "checking options.\n";

#checking for gzcat
my $zcat = "zcat";
my $gzcat_test = `which gzcat`;
if ($gzcat_test ne '') {
	$zcat = "gzcat";
}

#setting the default values.
my $r1 = ""; #file path to R1.
my $r2 = ""; #file path to R2.
my $index = ""; #library index sequence
my $out = ""; #output directory
my $mode = "perl"; #"fastqc_pair" or "perl"
my $threads = "2"; #number of threads given

my $raw_folder = "01_RAW";
my $trimmed_folder = "02_TRIMMED";
my $multi_folder = "03_MULTIQC";
my $integrated_folder = "04_INTEGRATED";

my $trimmings = "trimmed_sequences.fas";
my $m_threshold = 1; #shortest contaminant to consider
my $l_threshold = 15; #shortest length of sequence to keep after trimming
my $p_threshold = 0.99; #my prior
my $q_threshold = 30; #quality threshold

my $log = "TrimSeq4scRNAseq_log.txt";
my $outdir = "TrimSeq4scRNAseq_out";

#making the options into external arguments.
GetOptions (
	'r1=s' => \$r1,
	'r2=s' => \$r2,
	'index=s@' => \$index,
	'out=s' => \$out,
	'mode=s' => \$mode,
	'q_threshold=s' => \$q_threshold,
	'l_threshold=s' => \$l_threshold,
	'threads' => \$threads,
	);

#checking for required options.
if (!$r1) {
	die "USAGE: option --r1 is required.\n";
}
if (!$r2) {
	die "USAGE: option --r2 is required.\n";
}
if (!$index) {
	die "USAGE: option --index is required.\n";
}
if (!$out) {
	die "USAGE: option --out is required.\n";
}
my @indices = @{ $index };

#checking option quality
if (-e $out && -d $out) {
	$outdir = $out."/".$outdir."_Q".$q_threshold."L".$l_threshold;
	$log = $outdir."/".$log;
	system "rm -rf $outdir";
	system "mkdir $outdir";
	print " out directory checked.\n";
}
else {
	die "Error: Selected output directory $out does not exist.\n";
}

#grabbing versions of tools
my $fastqc_version = `fastqc --version | head -n 1 | cut -f2 -d' ' | perl -pe chomp`;
my $scythe_version = `scythe --version | head -n 1 | cut -f3 -d' ' | perl -pe chomp`;
my $sickle_version = `sickle se --version | head -n 1 | cut -f3 -d' ' | perl -pe chomp`;
my $multiqc_verison = `multiqc --version | cut -f3 -d' ' | perl -pe chomp`;

#open a log file
open (LOG, ">", "$log") or die "cannot open $log.\n";
print LOG "#####TrimSeq4scRNAseq.pl LOG#####\n";
print LOG "R1 file:\t$r1\n";
print LOG "R2 file:\t$r2\n";
print LOG "Index sequence:\t@indices\n";
print LOG "\n";
print LOG "Tool versions:\n";
print LOG "\tFASTQC $fastqc_version\n";
print LOG "\tscythe $scythe_version\n";
print LOG "\tsickle $sickle_version\n";
print LOG "\tmultiqc $multiqc_verison\n";
print LOG "\n";
##########



#####Step-01: Make input file folder#####
print "Step-01: Making raw input folder\n";
print LOG "Step-01: Making raw input folder\n";

$raw_folder = $outdir."/".$raw_folder;
system "mkdir $raw_folder";

my $r1_raw_file = `basename $r1 | perl -pe chomp`;
$r1_raw_file = $raw_folder."/".$r1_raw_file;
my $r2_raw_file = `basename $r2 | perl -pe chomp`;
$r2_raw_file = $raw_folder."/".$r2_raw_file;

system "cp $r1 $raw_folder";
system "cp $r2 $raw_folder";

if ($r1_raw_file !~ m/.gz$/) {
	print " compressing $r1_raw_file\n";
	system "gzip $r1_raw_file";
	$r1_raw_file = $r1_raw_file.".gz";
	print LOG "\tR1 file compressed.\n";
}
if ($r2_raw_file !~ m/.gz$/) {
	print " compressing $r2_raw_file\n";
	system "gzip $r2_raw_file";
	$r2_raw_file = $r2_raw_file.".gz";
	print LOG "\tR2 file compressed.\n";
}

print " running FASTQC on raw R1 and R2 files\n";
system "fastqc $r1_raw_file $r2_raw_file -o $raw_folder -f fastq -t $threads";
print LOG "\tcmd: fastqc $r1_raw_file $r2_raw_file -o $raw_folder -f fastq -t $threads\n";
##########



#####Step-02: run trimming of R2 file#####
print "Step-02: Running trimming on raw R2 file\n";
print LOG "Step-02: Running trimming on raw R2 file\n";
print LOG "Parameter settings:\n";
print LOG "\tShortest contaminant length:\t$m_threshold\n";
print LOG "\tShortest sequence length:\t$l_threshold\n";
print LOG "\tPrior probability:\t$p_threshold\n";
print LOG "\tQuality threshold:\t$q_threshold\n";

$trimmed_folder = $outdir."/".$trimmed_folder;
system "mkdir $trimmed_folder";

my $r2_at_file = `basename $r2_raw_file | perl -pe chomp`;
$r2_at_file = $trimmed_folder."/at.".$r2_at_file;

my $r2_atqt_file = `basename $r2_raw_file | perl -pe chomp`;
$r2_atqt_file = $trimmed_folder."/atqt.".$r2_atqt_file;

$trimmings = $trimmed_folder."/".$trimmings;

my $polyA = "AAAAAAAAAAAAAAA";
my $polyT = "TTTTTTTTTTTTTTT";

my $univ_adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
my $univ_adapter_rc = reverse $univ_adapter;
$univ_adapter_rc =~ tr/ATGC/TACG/;
my $se_pcr_primer1 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
my $se_pcr_primer1_rc = reverse $se_pcr_primer1;
$se_pcr_primer1_rc =~ tr/ATGC/TACG/;
my $rna_pcr_primer = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA";
my $rna_pcr_primer_rc = reverse $rna_pcr_primer;
$rna_pcr_primer_rc =~ tr/ATGC/TACG/;

my $nextera_transposase_1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
my $nextera_transposase_1_rc = reverse $nextera_transposase_1;
$nextera_transposase_1_rc =~ tr/ATGC/TACG/;
my $nextera_transposase_2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
my $nextera_transposase_2_rc = reverse $nextera_transposase_2;
$nextera_transposase_2_rc =~ tr/ATGC/TACG/;

my $clontech_univ_primer_long = "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT";
my $clontech_univ_primer_long_rc = reverse $clontech_univ_primer_long;
$clontech_univ_primer_long_rc =~ tr/ATGC/TACG/;
my $smart_cds_primer_ii_a = "AAGCAGTGGTATCAACGCAGAGTACT";
my $smart_cds_primer_ii_a_rc = reverse $smart_cds_primer_ii_a;
$smart_cds_primer_ii_a_rc =~ tr/ATGC/TACG/;
my $smarter_ii_a_oligo = "AAGCAGTGGTATCAACGCAGAGTAC";
my $smarter_ii_a_oligo_rc = reverse $smarter_ii_a_oligo;
$smarter_ii_a_oligo_rc =~ tr/ATGC/TACG/;

open (TRIMMINGS, ">", "$trimmings") or die "cannot open $trimmings\n";
print TRIMMINGS ">As\n";
print TRIMMINGS "$polyA\n";
print TRIMMINGS ">Ts\n";
print TRIMMINGS "$polyT\n";
print TRIMMINGS ">IlluminaUniversalAdapter\n";
print TRIMMINGS "$univ_adapter\n";
print TRIMMINGS ">IlluminaUniversalAdapter_rc\n";
print TRIMMINGS "$univ_adapter_rc\n";
print TRIMMINGS ">IlluminaSEPCRPrimer1\n";
print TRIMMINGS "$se_pcr_primer1\n";
print TRIMMINGS ">IlluminaSEPCRPrimer1_rc\n";
print TRIMMINGS "$se_pcr_primer1_rc\n";
print TRIMMINGS ">IlluminaRNAPCRPrimer\n";
print TRIMMINGS "$rna_pcr_primer\n";
print TRIMMINGS ">IlluminaRNAPCRPrimer_rc\n";
print TRIMMINGS "$rna_pcr_primer_rc\n";
print TRIMMINGS ">NexteraTransposaseRead1\n";
print TRIMMINGS "$nextera_transposase_1\n";
print TRIMMINGS ">NexteraTransposaseRead1_rc\n";
print TRIMMINGS "$nextera_transposase_1_rc\n";
print TRIMMINGS ">NexteraTransposaseRead2\n";
print TRIMMINGS "$nextera_transposase_2\n";
print TRIMMINGS ">NexteraTransposaseRead2_rc\n";
print TRIMMINGS "$nextera_transposase_2_rc\n";
foreach my $ind (@indices) {
	my $nextera_pcr_primer_i7 = "CAAGCAGAAGACGGCATACGAGAT".$ind."GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG";
	my $nextera_pcr_primer_i7_rc = reverse $nextera_pcr_primer_i7;
	$nextera_pcr_primer_i7_rc =~ tr/ATGC/TACG/;
	
	print TRIMMINGS ">NexTera_PCR_primer_i7_$ind\n";
	print TRIMMINGS "$nextera_pcr_primer_i7\n";
	print TRIMMINGS ">NexTera_PCR_primer_i7_rc_$ind\n";
	print TRIMMINGS "$nextera_pcr_primer_i7_rc\n";
}
print TRIMMINGS ">ClontechUniversalPrimerLong\n";
print TRIMMINGS "$clontech_univ_primer_long\n";
print TRIMMINGS ">ClontechUniversalPrimerLong_rc\n";
print TRIMMINGS "$clontech_univ_primer_long_rc\n";
print TRIMMINGS ">SMART_CDS_primer_IIA\n";
print TRIMMINGS "$smart_cds_primer_ii_a\n";
print TRIMMINGS ">SMART_CDS_primer_IIA_rc\n";
print TRIMMINGS "$smart_cds_primer_ii_a_rc\n";
print TRIMMINGS ">SMARTer_II_A_oligo\n";
print TRIMMINGS "$smarter_ii_a_oligo\n";
print TRIMMINGS ">SMARTer_II_A_oligo_rc\n";
print TRIMMINGS "$smarter_ii_a_oligo_rc\n";
close (TRIMMINGS);

#run scythe
print " running adapter trimming with scythe\n";
system "scythe -a $trimmings -q sanger -p $p_threshold -n $m_threshold -M $l_threshold --quiet $r2_raw_file | gzip >$r2_at_file";
print LOG "\tcmd: scythe -a $trimmings -q sanger -p $p_threshold -n $m_threshold -M $l_threshold --quiet $r2_raw_file | gzip >$r2_at_file\n";

#run sickle
print " running adapter trimming with sickle\n";
system "sickle se -f $r2_at_file -t sanger -q $q_threshold -l $l_threshold -g --quiet -o $r2_atqt_file";
print LOG "\tcmd: sickle se -f $r2_at_file -t sanger -q $q_threshold -l $l_threshold -g --quiet -o $r2_atqt_file\n";

#run fastqc
print " running fastqc on adapter and quality trimmed files\n";
system "fastqc $r2_at_file $r2_atqt_file -o $trimmed_folder -f fastq -t $threads";
print LOG "\tcmd: fastqc $r2_at_file $r2_atqt_file -o $trimmed_folder -f fastq -t $threads\n";
##########



#####Step-03: run multi fastq_pair#####
print "Step-03: Running multiqc for all fastqc files\n";
print LOG "Step-03: Running multiqc for all fastqc files\n";

$multi_folder = $outdir."/".$multi_folder;
system "mkdir $multi_folder";

system "mv $raw_folder/*fastqc* $multi_folder";
system "mv $trimmed_folder/*fastqc* $multi_folder";

#run multi qc
print " running multiqc on all fastqc outputs\n";
system "multiqc $multi_folder -o $multi_folder";
print LOG "\tcmd: multiqc $multi_folder -o $multi_folder\n";
##########



#####Step-04: pairing R2 with R1#####
print "Step-04: Running integration of R1 and R2\n";
print LOG "Step-04: Running integration of R1 and R2\n";

$integrated_folder = $outdir."/".$integrated_folder;
system "mkdir $integrated_folder";

my $r1_paired_file = `basename $r1_raw_file | perl -pe chomp`;
$r1_paired_file = $integrated_folder."/".$r1_paired_file;
my $r2_paired_file = `basename $r2_raw_file | perl -pe chomp`;
$r2_paired_file = $integrated_folder."/".$r2_paired_file;
my $r1_single_file = `basename $r1_raw_file | perl -pe chomp`;
$r1_single_file = $integrated_folder."/single.".$r1_single_file;

my $temp_r1 = $r1_raw_file.".temp";
my $temp_r2 = $r2_raw_file.".temp";
system "$zcat $r1_raw_file > $temp_r1";
system "$zcat $r2_atqt_file > $temp_r2";

my $temp_r1_paired = $temp_r1.".paired.fq";
my $temp_r2_paired = $temp_r2.".paired.fq";
my $temp_r1_single = $temp_r1.".single.fq";
my $temp_r2_single = $temp_r2.".single.fq";

print " integrating R1 and R2 files\n";
if ($mode =~ m/fastqc_pair/) {
	system "fastq_pair $temp_r1 $temp_r2";
	print LOG "\tcmd: fastq_pair $temp_r1 $temp_r2\n";
}
elsif ($mode =~ m/perl/) {
	my $linecount = 0;
	my $linecount_full = `wc -l $temp_r1 | rev | cut -f2 -d' ' | rev | perl -pe chomp`;
	
	open (TEMP1, "<", $temp_r1) or die "cannot open $temp_r1\n";
	open (TEMP2, "<", $temp_r2) or die "cannot open $temp_r2\n";
	open (PR1, ">", $temp_r1_paired) or die "cannot open $temp_r1_paired\n";
	open (SR1, ">", $temp_r1_single) or die "cannot open $temp_r1_single\n";
	while (my $h = <TEMP2>) {
		my $s = <TEMP2>;
		my $b = <TEMP2>;
		my $q = <TEMP2>;
		
		my $ref = (split (/ /, $h))[0];
		my $test = 0;
		
		while ($test == 0) {
			$linecount = $linecount + 4;
			if ($linecount > $linecount_full) {
				die "Error: $ref is not found in the untrimmed R1 file.\n";
			}
			
			my $header = <TEMP1>;
			my $seq = <TEMP1>;
			my $blank = <TEMP1>;
			my $qual = <TEMP1>;
			
			if ($header =~ m/$ref/) {
				$test = 1;
				print PR1 "$header";
				print PR1 "$seq";
				print PR1 "$blank";
				print PR1 "$qual";
			}
			else {
				print SR1 "$header";
				print SR1 "$seq";
				print SR1 "$blank";
				print SR1 "$qual";
			}
		}
	}
	close (TEMP1);
	close (TEMP2);
	close (PR1);
	close (SR1);
	
	my $tail = $linecount_full - $linecount;
	system "tail -n $tail $temp_r1 >> $r1_single_file";
}

system "gzip -c $temp_r2 > $r2_paired_file";
system "gzip -c $temp_r1_paired > $r1_paired_file";
system "gzip -c $temp_r1_single > $r1_single_file";

unlink ($temp_r1_paired);
unlink ($temp_r2_paired);
unlink ($temp_r1_single);
unlink ($temp_r2_single);

unlink ($temp_r1);
unlink ($temp_r2);
##########



#####Step-05: completing process#####
print "Step-05: Completing process\n";
print LOG "Step-05: Completing process\n";
print LOG "##########\n";
close (LOG);
##########
