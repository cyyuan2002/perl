#!/usr/bin/perl
use strict;

use strict;
use Getopt::Long;

my $version="1.0 Alvin Chen 2011-02-18";
my %opts;
##Read parameters
GetOptions(\%opts,"s=s","i=s","o:s","m:i","help");
if((!defined $opts{s})||(!defined $opts{i})){
    &Usage();
}

my $parSequenceFile=$opts{s};
my $parTranslateFile=$opts{i};
my $parOutputFile=(defined $opts{o}) ? $opts{o} : "$parTranslateFile.lgt";
my $parMinLength=(defined $opts{c}) ? $opts{c} : 0;

my $iserror=0;

&Filecheck($parSequenceFile,0);
&Filecheck($parTranslateFile,0);

if($iserror==1){
    exit(1);
}

my $seqfilein;
my $transfilein;
my $outfile;
open ($seqfilein, $parSequenceFile);
my @seqNames;
while(<$seqfilein>){
    chomp();
    if(/^>/){
	push (@seqNames,$_);
    }
}
close $seqfilein;

open ($outfile,">$parOutputFile");
open($transfilein,$parTranslateFile);
my $lastSeqN;
my $lastSeq;
my $lastSeqLength;
my $lastSeqSite;
my $SeqSite;
my $seqID;
my $seqcount=-1;
while(<$transfilein>){
    chomp();
    if(/^>(\S+)/){
	my @seqN=split(/\_/,$1);
	my $seqNs=$seqN[0];
	my $seqinfo=$_;
	if($seqinfo=~/(\[\d+\s-\s\d+\])/){
	    $SeqSite=$1;
	}
	if($seqNs ne $lastSeqN){
	    if($lastSeqN ne ""){
		$seqID=$seqNames[$seqcount];

		if(!($seqID=~/$lastSeqN/)){
		    print stderr "Error!";
		    exit(1);
		}
		if($parMinLength==0){
		    print $outfile "$seqID $lastSeqSite\n$lastSeq\n";
		}
		else{
		    if($lastSeqLength>$parMinLength){
			print $outfile "$seqID $lastSeqSite\n$lastSeq\n";
		    }
		}
	    }
	    $seqcount++;
	    $lastSeqN=$seqNs;
	    $lastSeq="";
	    $lastSeqLength="";
	}
    }
    else{
	if($lastSeq eq ""){
	    $lastSeq=$_;
	    $lastSeqSite=$SeqSite;
	    $lastSeqLength=length($_);
	}
	else{
	    if($lastSeqLength<length($_)){
		$lastSeq=$_;
		$lastSeqSite=$SeqSite;
		$lastSeqLength=length($_);
	    }

	}
    }
}

$seqID=$seqNames[$seqcount];
if($parMinLength==0){
    print $outfile "$seqID $lastSeqSite\n$lastSeq\n";
}
else{
    if($lastSeqLength>$parMinLength){
	print $outfile "$seqID $lastSeqSite\n$lastSeq\n";
    }
}


sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "File $Filename already exists\n";
	    $iserror=1;
	}
    }
}

sub Usage(){
  print << "    Usage";

	This program is used to filter the result of EMBOSS getorf program. The result is the longest protein sequences.
	The program is optimised for trans-abyss assembled sequences.

	Usage:  $0 (version $version)

	<options>
		-s     Sequence file
		-i     EMBOSS getorf output file
		-o     Output file (default: getorf.lgt)
		-m     Minimal protein length (default: off)

    Usage

	exit(0);
};
