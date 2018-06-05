#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#! /usr/bin/perl

## General info:
# script for analysing SAM files containing aligned results from 4C experiments.
# script will cycle through SAM results, will match read1 to read2, count read2 for each read 

# input files: 
# 		1. read 1 alignent result (containing all digestion products, to jump over non-cut fragments)
#		2. read 2 alignment result (only the primary digestion product)
# Output: compiled table of ligation events. 
# 		1. unique table of ligation sites (unique R1)
# 		2. unique table of ligation events (unique R1-R2 pairs)
# 		3. Summary table for extraction

# run example: 
#			perl FC04_analyseSAMmapq.pl <r1file.sam.gz> <r2file.sam.gz> \
# 								<out.ligsites> <out.ligevents> <out.summary>


use warnings;
use strict;

########################################################################
## subroutines

#################################################
## Reverse complement

sub revcom {

my ($DNA) = @_;
my $revDNA = reverse ($DNA);
$revDNA =~ tr/[]ACGTacgt/][TGCAtgca/;
return $revDNA

}

#################################################
## calculate binary number

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0{21}//;   # otherwise you'll get leading zeros
    return $str;
}

#################################################
## add new hash reference to site hashs

sub add_newsite {
    my($hashrefname, $hashname, $chromosome, $position, $strand, $SeqExample, $mapq, $AS, $diff ) = @_;
    $hashrefname->{$hashname} = {
        chromosome  => $chromosome,
        position => $position,
	strand => $strand,
	SeqExample => $SeqExample,
	ligevent => 0,
	ambligevent => 0, 
	mapq => $mapq, 
	AS => $AS, 
	diff => $diff,
	ligread => 1
    };
}
#################################################
## add new ligation event to hash. 

sub add_newevent {
    my($hashrefname, $hashname, $chromosome, $position, $strand, 
    $shears, $SeqExample, $Seq2Example, $mapq, $AS, $diff ) = @_;
    $hashrefname->{$hashname} = {
        chromosome  => $chromosome,
        position => $position,
	strand => $strand,
	shearsite => $shears, 
	SeqExample => $SeqExample,
	Seq2Example => $Seq2Example,
	mapq => $mapq,
	AS => $AS, 
	diff => $diff,
	ligread => 1
    };
}

#
sub add_shearsite {
    my($hashrefname, $hashname, $shearsite ) = @_;
    $hashrefname->{$hashname} = {
	$shearsite => 1
    };
}


################################################
## parse cigar sequence and calculate actual start position in reverse stranded seq. 

sub parsecigar {
	### cigar parsing modified from http://www.bioinformatics.babraham.ac.uk/projects/bismark/deduplicate_bismark ()

	my( $cigarstring, $position ) = @_;
	
	$position -= 1;

	
	# for InDel free matches we can simply use the M number in the CIGAR string
	if ($cigarstring =~ /^(\d+)M$/){ # linear match
	$position += $1;
	}

	else{
	# parsing CIGAR string
		my @len = split (/\D+/,$cigarstring); # storing the length per operation
		my @ops = split (/\d+/,$cigarstring); # storing the operation
		shift @ops; # remove the empty first element
		die "CIGAR string contained a non-matching number of lengths and operations\n" unless (scalar @len == scalar @ops);
#				warn "CIGAR string; $cigarstring\n";
		### determining end position of a read
		foreach my $index(0..$#len){
			
			if ($ops[$index] eq 'M'){  # standard matching bases
				$position += $len[$index];
			#	warn "Operation is 'M', adding $len[$index] bp\n";
			}
			
			elsif($ops[$index] eq 'I'){ # insertions do not affect the end position
			#	warn "Operation is 'I', next\n";
			}
			
			elsif($ops[$index] eq 'D'){ # deletions do affect the end position
			#	warn "Operation is 'D',adding $len[$index] bp\n";
				$position += $len[$index];
			}
			
			else{
			die "Found CIGAR operations other than M, I or D: '$ops[$index]'. Not allowed at the moment\n";
			}
		}
	}
	return $position;
}

########################################################################
## script

## get inputfile, outputfile names from command line input
my ($inputfilename, $inputfilename2, $output1filename, $output2filename, $summaryfilename)=@ARGV;


## declare variables
my @samline = '';
my $startbit = '';
my $location = '';
my $rname = '';
my $rpos = 0;
my $readlength = 0;
my $readsite = 0;
my $ligsite = '';
my $ligevent = '';
my $shearsite = '';
my $strand = '';
my $read1seq = '';
my $read2seq = '';
my $mapq1 = '';
my $flag = '';
my $AS1 = '';
my $diff1 = '';
my $XS1 = '';
my @flagfull = '';
my %shearhash = ();
my %shearseq = ();
my %shearhash2 = ();
my %ligsitehash = ();
my %ligeventhash = ();
my $cigar = '';


## declare counters
my $counter = 0;
my $fragcount = 0;
my $fragnouse = 0;
my $frag = 0;
my $frag2do = 0;
my $nomatch = 0;
my $ligevents = 0;
my $ambig = 0;

## test input:
scalar( @ARGV ) == 5 || die( "Analysis of paired SAM files requires input of 5 file names\n" );


######################################
######################################
## loop across SAM file starts here ## 
######################################
######################################


##############
# read read2 #
##############

# read2
unless ( open (IN2, "gunzip -c $inputfilename2 | " ) ){
	print "fileopen failed, please check $inputfilename2!\n";
	exit;
}

# report to user
print "Reading Read2 file, $inputfilename2\n";

while (<IN2>){
	
	# report to user
	if ($.%100000 == 0) {
	  print ".";
	}
	
	#remove file header
	if (/(^\@.*)/){
	next;
	}
	chomp;
	
	# record read details
	@samline=split(" ", $_);
	$location = $samline[0];
	$flag = $samline[1];
	@flagfull=split( "", dec2bin($flag) );
	$rname = $samline[2];
	$rpos = $samline[3];
	$cigar = $samline[5];
	$readlength = length( $samline[9] );
	if ( $flagfull[6] ) {
		$read2seq = revcom($samline[9]);
		$readsite=parsecigar( $cigar, $rpos ); 
	}
	else{		
		$read2seq = $samline[9];
		$readsite=$rpos;
	}
	$shearsite = $rname."_".$readsite."_".$flag;
	if ( $samline[2] eq "*" || $samline[4]<10 ) {
	$shearsite = "***";
	}
	# record shear site position of the cluster
	$shearhash{$location} = $shearsite;
	
	# record an example sequence for each shear site.
	# NOTE: the sequence sample could theoretically come from > lig site.
	# (shouldnt matter if is )
	if ( !exists($shearseq{$shearsite}) ){
		$shearseq{$shearsite} = $read2seq;
	}else{
		if ( length $read2seq > length $shearseq{ $shearsite } ) {
		    $shearseq{$shearsite} = $read2seq;
		}
	}


}
close IN2;


###################
# reset variables #
###################
@samline='';
@flagfull='';
$location = '';
$flag = '';
$rname = '';
$rpos = 0;
$readlength = 0; 
$readsite='';
$shearsite='';
$cigar='';

##############
# read read1 #
##############

# read1
unless ( open (IN1, "gunzip -c $inputfilename |" ) ){
	print "fileopen failed, please check $inputfilename!\n";
	exit;
}

# report to user
print "\nReading Read1 file, $inputfilename\n";

while (<IN1>){

	# report to user
	if ($.%100000 == 0) {
	  print ".";
	}

######### reset variable
	$AS1 = $XS1 = $diff1 = $mapq1 = "";
	
######### remove file header
	if (/(^\@.*)/){
	next;
	}
	chomp;

	$fragcount++;
	
######### record read details
	@samline=split(" ", $_);
	
	$startbit = $samline[0];
	$location = $startbit;
	$location=~s/:\d$//g;
	$frag=substr( $startbit, (length($startbit))-1 );

	$AS1 = $samline[11];
	$AS1=~s/AS:i://g;

	$XS1 = $samline[12];

	$flag = $samline[1];
	@flagfull=split( "", dec2bin($flag) );
	$rname = $samline[2];
	$rpos = $samline[3];
	$mapq1 = $samline[4];

	$readlength = length( $samline[9] );

	$cigar = $samline[5];
	if ( $cigar =~ /(^\d+)M/) {
	}

	

	

######## which fragment to do? 
	if ( $frag == 1 && $rname eq "HTLV-1" && $rpos == 7243  ) {
		$frag2do=$location.":2";
		$fragnouse++;
		next;
	}
	if ( $frag == 1 && !( $rname eq "HTLV-1" && $rpos == 7243 ) ) {
		$frag2do=$startbit;
	}
	if ( $frag > 1 && $frag2do ne $startbit ){
		$frag2do=0;
		$fragnouse++;
		next;
	}

	$counter++;

	if ( $samline[2] eq "*" || $samline[4]<10 || $AS1 <= -10 ) {
		$nomatch++;
		next;
	}

	
######### if sequence is mapped to forward strand, keep as is. if on reverse strand -
######### reverse complement sequence. 
	if ( $flagfull[6] ) {
		$read1seq = revcom($samline[9]);
		$strand = "-";
#		$readsite=$rpos+$readlength-1;
		$readsite=parsecigar( $cigar, $rpos ); 
	}
	else{		
		$read1seq = $samline[9];
		$strand = "+";
		$readsite=$rpos;
	}

	$ligsite = $rname."_".$readsite."_".$flag;
	$shearsite = $shearhash{$location};
	$ligevent=$ligsite.":".$shearsite;

	if ($XS1 =~ /XS/)
	{
		$XS1=~s/XS:i://g;
		$diff1=$AS1-$XS1;
	} else {
		$diff1 = 999;
	}


	if ( !exists($shearhash{$location}) ){
	print "shear site information missing: location is $location\t";
	print "shearsite is $shearsite\n";
	}
	
	
######### if this ligation event is observed for the first time, add to the ligeventhash
######### if has been observed, check if can add longer sequence example. 

	if ( !exists($ligeventhash{$ligevent}) ){
		add_newevent \%ligeventhash, $ligevent, $rname,
				$readsite, $strand, $shearsite,
				$read1seq, $shearseq{$shearsite}, $mapq1, $AS1, $diff1;
	}else{
		$ligeventhash{$ligevent}{ligread}++;
		if ( length($ligeventhash{$ligevent}{SeqExample})<length($read1seq) ) {
			$ligeventhash{$ligevent}{SeqExample}=$read1seq
		}
		if ( $mapq1 < $ligeventhash{$ligevent}{'mapq'} ) {
			$ligeventhash{$ligevent}{'mapq'} = $mapq1;
		}
		if ( $AS1 < $ligeventhash{$ligevent}{'AS'} ) {
			$ligeventhash{$ligevent}{'AS'} = $AS1;
		}
		if ( $diff1 < $ligeventhash{$ligevent}{'diff'} ) {
				$ligeventhash{$ligevent}{'diff'} = $diff1;
		}
	}

######### if this ligation site are observed for the first time, add to the ligsitehash
######### if has been observed, check if can add longer sequence example. 

	if ( !exists($ligsitehash{$ligsite}) ){
		add_newsite \%ligsitehash, $ligsite, $rname, $readsite, 
		$strand, $read1seq, $mapq1, $AS1, $diff1;
	}else{
		$ligsitehash{$ligsite}{ligread}++;
		if ( length($ligsitehash{$ligsite}{SeqExample})<length($read1seq) ) {
			$ligsitehash{$ligsite}{SeqExample}=$read1seq
		}
		if ( $mapq1 < $ligsitehash{$ligsite}{'mapq'} ) {
			$ligsitehash{$ligsite}{'mapq'} = $mapq1;
		}
		if ( $AS1 < $ligsitehash{$ligsite}{'AS'} ) {
			$ligsitehash{$ligsite}{'AS'} = $AS1;
		}
		if ( $diff1 < $ligsitehash{$ligsite}{'diff'} ) {
			$ligsitehash{$ligsite}{'diff'} = $diff1;
		}
	}

######### if this shear site does not exit or unmapped, list as ambigous lig event
######### if this shear site exists: 
######### if this shear site is observed for the first time for this , add to the ligsitehash and update the ligation site event count 
######### if has been observed, add to the count of the shear site. 

	if ( !exists($shearhash{$location}) ) {
		$ligsitehash{$ligsite}{ambligevent}++;
		$ambig++;
	}elsif ( $shearsite eq "***" ) {
		$ligsitehash{$ligsite}{ambligevent}++;
		$ambig++;
	}else{
		if ( !exists($shearhash2{$ligsite}) ){
			add_shearsite \%shearhash2, $ligsite, $shearsite;
			$ligsitehash{$ligsite}{ligevent}++;
			$ligevents++;
		}elsif( exists($shearhash2{$ligsite}{$shearsite}) ) {
			$shearhash2{$ligsite}{$shearsite}++;
		}else{
			$shearhash2{$ligsite}{$shearsite} = 1; 
			$ligsitehash{$ligsite}{ligevent}++;
			$ligevents++;
		}
	}
}
close IN1;

#############################################
#############################################
## output to user - list of ligation sites ##
#############################################
#############################################

print "\nWriting Output files\n";

unless ( open (OUT1, ">$output1filename" ) ){
	print "fileopen failed, please check $output1filename!\n";
	exit;
}

print OUT1 "SiteID\tChr\tPosition\tStrand\tSeqExample\tLigEvents\tAmbLigReads\tmapq\tAS\tdiff\tReadCount\n";
my ($l)='';
for (sort keys %ligsitehash){
        $l=$_;
        print OUT1 $_,
	"\t", $ligsitehash{$_}{chromosome},
	"\t", $ligsitehash{$_}{position},
	"\t", $ligsitehash{$_}{strand},
	"\t", $ligsitehash{$_}{SeqExample},
	"\t", $ligsitehash{$_}{ligevent},
	"\t", $ligsitehash{$_}{ambligevent},
	"\t", $ligsitehash{$_}{mapq},
	"\t", $ligsitehash{$_}{AS},
	"\t", $ligsitehash{$_}{diff},
	"\t", $ligsitehash{$_}{ligread},
	"\n";
}
close OUT1;

##############################################
##############################################
## output to user - list of ligation events ##
##############################################
##############################################
unless ( open (OUT2, ">$output2filename" ) ){
	print "fileopen failed, please check $output2filename!\n";
	exit;
}

print OUT2 "EventID\tChr\tPosition\tStrand\tSeqExample\tSeq2Example\tmapq\tAS\tdiff\tReadCount\n";

for (sort keys %ligeventhash){
        print OUT2 $_,
	"\t", $ligeventhash{$_}{chromosome},
	"\t", $ligeventhash{$_}{position},
	"\t", $ligeventhash{$_}{strand},
	"\t", $ligeventhash{$_}{SeqExample},
	"\t", $ligeventhash{$_}{Seq2Example},
	"\t", $ligeventhash{$_}{mapq},
	"\t", $ligeventhash{$_}{AS},
	"\t", $ligeventhash{$_}{diff},
	"\t", $ligeventhash{$_}{ligread},
	"\n";
}


###################################
###################################
## output to user - summary file ##
###################################
###################################

unless ( open (OUT3, ">$summaryfilename" ) ){
	print "fileopen failed, please check $summaryfilename!\n";
	exit;
}

	print OUT3 "Total read1 fragments: $fragcount\n";
	print OUT3 "Total read1 fragments unused (non-cut or frag>1 where frag1 was not non-cut): $fragnouse\n";	
	print OUT3 "Total reads: $counter\n";
	print OUT3 "Total not uniquely aligned: $nomatch\n";
	print OUT3 "Total reads uniquely aligned:", ($counter-$nomatch), "\n";
	print OUT3 "Total ligation events across all ligation sites: $ligevents\n";
	print OUT3 "Total ambiguous ligation reads across all ligation sites: $ambig\n";
	print OUT3 "Total unique ligation sites: ",scalar keys( %ligsitehash ), "\n";

	close OUT3;

exit;

