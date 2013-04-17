#!/usr/bin/perl

use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","s:s","l:s","help");

if(!defined $opts{i}){
    &Usage();
}
my $inputfile=$opts{i};
my $outputfile=(defined $opts{o}) ? $opts{o} :"$inputfile.flt";
my $scorecut=(defined $opts{s}) ? $opts{s} : 20;
my $lengthcut=(defined $opts{l}) ? $opts{l} : 17;

open(file,$inputfile) || die "Can't open file: $inputfile\n";


#print stderr "Reading File... $inputfile\n";
my @reads=<file>;
close file;
#print stderr "Finished Reading...Processing\n";

open (outfile,">$outputfile");


for(my $i=0;$i<@reads;$i=$i+4){
    if(length($reads[$i+1])<36){
        next;
    }
    my $sageseq=substr($reads[$i+1],0,$lengthcut);
    if($sageseq=~/N/){
        next;
    }
    my $scores=substr($reads[$i+3],0,$lengthcut);
    my @letters=split("",$scores);
    my $ispass=1;
    for(my $j=0;$j<17;$j++){
        my $score=ord($letters[$j])-64;
        if($score<$scorecut){
            $ispass=0;
            last;
        }
    }
    if($ispass==1){
        print outfile "$reads[$i]$sageseq\n$reads[$i+2]$scores\n";
    }
}
close outfile;

