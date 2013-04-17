#!/usr/bin/perl
use strict;
use Getopt::Long;

my %opts;
##Read parameters
GetOptions(\%opts,"i=s","o:s","f:i","m:i","p:i","l:i","a:i","e:s","s:s","help");
if(!defined $opts{i}){
    &Usage();
}

my $parInputFile=$opts{i};
my $parOutputFile=(defined $opts{o}) ? $opts{o} : "$parInputFile.rgn";
my $parOverlap=(defined $opts{p}) ? $opts{p} : "NULL";
my $parMatchLength=(defined $opts{l}) ? $opts{l} : "NULL";
my $parMatchPercent=(defined $opts{m}) ? $opts{m} : "NULL";
my $parIdentity=(defined $opts{a}) ? $opts{a} : "NULL";
my $parEvalue=(defined $opts{e}) ? $opts{e} : "NULL";
my $parScore=(defined $opts{s}) ? $opts{s} : "NULL";
my $parFlankLength=(defined $opts{f}) ? $opts{f} : "NULL";

my $iserror=0;
&Filecheck($parInputFile,0);
if($iserror==1){
    exit(1);
}

open(fileIN,$parInputFile);
my $lastID;
my $lastScore;
open (fileOut,">$parOutputFile");

my $seqcount=0;

while(<fileIN>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastID ne $lines[0]){
	$lastID = $lines[0];
	$lastScore = $lines[7];
	$seqcount=0;
    }
    my $left=$lines[4];
    my $right=$lines[5];
    my $dir="+";
    if($right<$left){
	my $mid=$left;
	$left=$right;
	$right=$mid;
	$dir="-";
    }
    my $size=$lines[6];
    if($lines[7]==$lastScore){
	if($parFlankLength ne "NULL"){
	    $left=(($lines[4]-$parFlankLength)<1) ? 1 : $lines[4]-$parFlankLength;
	    $right=(($lines[5]+$parFlankLength)>$size) ? $size : $lines[5]+$parFlankLength;
	}
	if($parOverlap ne "NULL"){
	    my @overlap=split(/\//,$lines[9]);
	    my $overpercent=$overlap[0]/$overlap[1];
	    if($overpercent<$parOverlap){
		next;
	    }
	}
	if($parMatchLength ne "NULL"){
	    my $matchlength=$lines[3]-$lines[2];
	    if($matchlength<$parMatchLength){
		next;
	    }
	}
	if($parMatchPercent ne "NULL"){
	    my $matchpercent=(($lines[3]-$lines[2]+1)/$lines[1])*100;
	    if($matchpercent<$parMatchPercent){
		next;
	    }
	}
	if($parIdentity ne "NULL" && $parIdentity>$lines[10]){
	    next;
	}
	if($parEvalue ne "NULL" && $parEvalue < $lines[8]){
	    next;
	}
	if($parScore ne "NULL" && $parScore>$lines[7]){
	    next;
	}
	print fileOut "$lines[0]\@$seqcount\t$lines[11]\t$dir\t$left\t$right\n";
	$seqcount++;
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

	Usage:  $0

	<options>
		-i     Input EblastN.pl output file
		-o     Output file, (default, Input_File.rgn)
		-f     Flanking length
		-p     Percentage of overlap
		-l     Minimal length of blast match
		-m     Minimal Match length percentage, (default 0.9)
		-a     Minimal identity
		-e     Maximum evalue
		-s     Minimal score


    Usage

	exit(0);
};
