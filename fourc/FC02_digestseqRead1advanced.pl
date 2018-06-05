#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#! /usr/bin/perl

## General info:
# script for digesting sequences in fastq file:
# read in each line in sets of 4
# if line 1: record information in scalar
# if line 2: sequence should start with CATG. if not move to next sequence.
# only read 1: require seqeunce start with CATG. 
# split according to CATG, record in array

# if a fragment is <35 try to add it back to the following fragment

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

my $digestedfilename = $stemname."-NlaIII.fq";
my $digestedfilename1 = $stemname."-first-NlaIII.fq";
my $pattern = "CATG";

## declare variables
my $line4 = '';
my $frag = 0; 
my $fragcount = 0; 
my $frags = 0; 
my $len1 = 0;
my $start = 0;
my $seq='';
my $startbit = '';
my $endbit = '';
my $header = '';
my $filter = '';
my $passfilter = 0;
my @digseq = '';
my @digseqover20 = ();
my $fragseq = '';
my $shortseq = '';

## declare variable for counting
my $counter = 0;
my $have4base = 0;
my $havepat = 0;
my $nopat = 0;
my $totalfrags = 0;
my $shortfrag = 0;
my $totalfragsfirst = 0;
my $shortfragfirst = 0;

## load sequence and output files.
## OUT1 - all fastq reads. OUT3 - only primary fragment

unless ( open (IN1, $inputfilename ) ){
	print "fileopen failed, please check $inputfilename!\n";
	exit;
}

unless ( open (OUT1, ">$digestedfilename" ) ){
	print "fileopen failed, please check $digestedfilename!\n";
	exit;
}

unless ( open (OUT3, ">$digestedfilename1" ) ){
	print "fileopen failed, please check $digestedfilename1!\n";
	exit;
}


## update to user
print scalar localtime(), " Digesting fastq file based on $pattern \n";

## loop across fastq file starts here
while (<IN1>){
	chomp;
	if (/(^$rundetails:.*)\s+(\w+:\w+:\w:\w+):*(\w*).*/){
		$startbit = $1;
		$endbit = $2;
		$filter = $3;
		print "$startbit..$endbit..$filter\n";
	next;
	}
 
 ## sequence retained only if starts with restriction site
 ## if so, digested such that only contains one restriction site
 ## for each fragment
  	if ($.%4 == 2){  		
		$counter++;    
		$len1 = 0;
		$passfilter=0;
		@digseq='';
		$frags=0;
		if (length ($filter) > 0 && $filter eq "CGGG") {
			$have4base++;
			$passfilter=1;
		}
		if (length ($filter) == 0) {
			$have4base="NA";
			$passfilter=1;
		}
		if ($passfilter) {
			if (/(^$pattern)/){
				$seq=$_;
				$seq =~ s/\s//g;                                                                   
				@digseq = split( $pattern, substr( $seq, 4, length($seq) ) ); 
				$frags = scalar( @digseq );
				$havepat++;
				next;
			}else{	
				$frags=0;			
				$nopat++;
				next;	
			}
		}else{
			next;
		}
	}

	if(/^\+$/){
		next;
	}

	if ($.%4 == 0){
		$line4 = $_;
	}
    
## setup output
	$frag = 0;
	$start = 0;
	$fragcount = 0; 
	@digseqover20 = ();
	
## process digseq to remerge fragments if < 20 base
	while ($fragcount < $frags) {
		push( @digseqover20, $digseq[$fragcount] );
		if ( scalar @digseqover20 > 1 ) {
			if ( length( $digseqover20[-2] ) < 20 ) {
				pop @digseqover20;
				$shortseq = pop @digseqover20;
				push( @digseqover20, $shortseq."CATG".$digseq[$fragcount] );
			} else {
			#	print "ok " ;
				}
		}
		$fragcount++;
	}
	$frags = scalar @digseqover20;

	$totalfrags = $totalfrags+$frags;

	
	while ($frag < $frags){
		$fragseq=$pattern.$digseqover20[$frag];
		if ( $frag < $frags-1 ) {
			$fragseq=$fragseq.$pattern;
		}

		$len1 = length( $fragseq );
		print OUT1 $startbit.":".($frag+1)." ".$endbit.":".$filter, "\n";
		print OUT1 $fragseq, "\n";
		print OUT1 "+\n";
		print OUT1 substr($line4, $start, $len1), "\n";
## OUT3 - only primary sequence (first object in @digseqover20 array)
		if ($frag == 0) {
			$totalfragsfirst++;
			print OUT3 $startbit.":".($frag+1)." ".$endbit.":".$filter, "\n";
			print OUT3 $fragseq, "\n";
			print OUT3 "+\n";
			print OUT3 substr($line4, $start, $len1), "\n";	
			if ($len1<35) {
				$shortfragfirst++;
			}
		}
		$start = $start+$len1-4;
		$frag++;
		if ($len1<35) {
			$shortfrag++;
		}
	}
## update to user to give feedback for very large files
	if ( $counter % 100000 == 0) {
		print ".";
	}
}   

close IN1; 
close OUT1; 
close OUT3; 

## open summary stats file and report seqeunces counts.
unless (-e $outputfilename) {
	$header="FileName\tTotClusters\tTotCGGG\tTotHave".$pattern."\tTotNo".$pattern."\tTotFrag"."\tTotShortFrag"."\tPctShortFrag"."\tPriFrag"."\tPriShortFrag"."\tPctPriShort\n";
} 
unless ( open (OUT2, ">>$outputfilename" ) ){
	print "fileopen failed, please check $outputfilename!\n";
	exit;
}

print OUT2 $header;
print OUT2 $inputfilename."\t".$counter."\t".$have4base."\t".$havepat."\t".$nopat."\t".$totalfrags."\t".$shortfrag."\t".100*(sprintf "%.4f", $shortfrag/$totalfrags )."\t".$totalfragsfirst."\t".$shortfragfirst."\t".100*(sprintf "%.4f", $shortfragfirst/$totalfragsfirst )."\n";


## update to user
print scalar localtime(), " done!\n";
 
exit;