#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

## General info:
# script for digesting sequences in fastq file:
# read in each line in sets of 4
# if line 1: record information in scalar
# if line 2: sequence should start with CATG. if not move to next sequence.
# only read 2: does not require seqeunce start with CATG and only keep
# first digest. if start with CATG it will be removed, so that is added back. 
# split according to CATG, record in array


use warnings;
use strict;

## get inputfile, outputfile names from command line input

my ($inputfilename, $outputfilename, $rundetails)=@ARGV;
my $stemname = $inputfilename;
$stemname =~ s/\.fq|\.fastq//g;

## check input
if( scalar @ARGV != 3 )	{
	die("enter input and output file names only and run name\n");
}

my $digestedfilename1 = $stemname."-first-NlaIII.fq";
my $pattern = "CATG";

## declare variables
my $line4 = '';
my $frags = 0; 
my $len1 = 0;
my $start = 0;
my $seq='';
my $startbit = '';
my $endbit = '';
my @digseq = '';
my $seqend = '';
my $firstfrag = '';
my $header = '';

## declare variable for counting
my $counter = 0;
my $totalfrags = 0;
my $totalfragsfirst = 0;
my $shortfragfirst = 0;

## load sequence and output files.
## OUT3 - only primary fragment

unless ( open (IN1, $inputfilename ) ){
	print "fileopen failed, please check $inputfilename!\n";
	exit;
}

unless ( open (OUT3, ">$digestedfilename1" ) ){
	print "fileopen failed, please check $digestedfilename1!\n";
	exit;
}

## update to user
print scalar localtime(), " Digesting fastq file based on $pattern \n";

print $rundetails, "\n";

## loop across fastq file starts here
while (<IN1>){
	chomp;
	if (/(^$rundetails:.*)\s+(\w+:\w+:\w:\w+).*/){
		$startbit = $1;
		$endbit = $2;
	next;
	}
 
 ## sequence retained only if starts with restriction site
 ## if so, digested such that only contains one restriction site
 ## for each fragment
  	if ($.%4 == 2){  		
		$counter++;    
		$len1 = 0;
		$seq=$_;
		$seq =~ s/\s//g;                                                                   
		@digseq = split( $pattern, $seq );
		$frags = scalar( @digseq );
		$seqend = substr( $seq, length($seq)-4 );
		if ( length $digseq[0] == 0 ) {
			shift @digseq;
			$digseq[0] = $pattern.$digseq[0];
		}

		if ( scalar @digseq == 1 && $seqend eq $pattern ) {
			$digseq[0] = $digseq[0].$pattern;
		}

		if ( scalar @digseq > 1 ) {
			$digseq[0] = $digseq[0].$pattern;
		}
		
		$firstfrag = $digseq[0];
	
		next;
	}

	if(/^\+$/){
		next;
	}

	if ($.%4 == 0){
		$line4 = $_;
	}
    
## setup output
	$totalfrags = $totalfrags+$frags;
		
## OUT3 - only final sequence (first object in @digseq array)
	$totalfragsfirst++;
	$len1 = length( $firstfrag );
	print OUT3 $startbit." "."1".":".$endbit, "\n";
	print OUT3 $firstfrag, "\n";
	print OUT3 "+\n";
	print OUT3 substr($line4, 0, $len1), "\n";	
	if ($len1<35) {
		$shortfragfirst++;
	}

	
## update to user to give feedback for very large files
	if ( $counter % 100000 == 0) {
		print ".";
	}
}   

close IN1; 
close OUT3; 

## open summary stats file and report seqeunces counts. 
unless (-e $outputfilename) {
	$header="FileName\tTotClusters\tTotFrag"."\tPriFrag"."\tPriShortFrag"."\tPctPriShort\n";
} 
unless ( open (OUT2, ">>$outputfilename" ) ){
	print "fileopen failed, please check $outputfilename!\n";
	exit;
}

print OUT2 $header;
print OUT2 $inputfilename."\t".$counter."\t".$totalfrags."\t".$totalfragsfirst."\t".$shortfragfirst."\t".100*(sprintf "%.4f", $shortfragfirst/$totalfragsfirst )."\n";

#print OUT2 "Total clusters: $counter\n";
#print OUT2 "Total sequences after digestion by $pattern: $totalfrags\n";
#print OUT2 "Total primary target sequences after digestion by $pattern:	$totalfragsfirst\n";
#print OUT2 "\tTotal primary target sequences under 35 bases: $shortfragfirst (",
#	100*(sprintf "%.4f", $shortfragfirst/$totalfragsfirst ),"%)\n";

close OUT2;

## update to user
print scalar localtime(), " done!\n";
 
exit;