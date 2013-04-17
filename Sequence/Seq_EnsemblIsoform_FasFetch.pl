#!/usr/bin/perl
##this program is used to get large number of fasta sequence in one time.

use strict;

my ($fasfile,$seqname,$outfile)=@ARGV;

if(@ARGV<3){
    print "Usage:$0 fastafile seqidfile outfile\n";
    exit;
}

my %seqnames;

open(seqN,$seqname) || die "Can't open file $seqname\n";
while(<seqN>){
    chomp();
    $seqnames{$_}=1;
}
close seqN;

my $isoutput=0;

open (outfile,">$outfile");

open(fasfile,$fasfile) || die "Can't open file $fasfile\n";
while(<fasfile>){
    chomp();
    if(/^>/){
        $seqname=$_;
        $seqname=~s/>//g;
        if(exists($seqnames{$seqname})){
            $isoutput=1;
            print outfile "$_\n";
        }
        else{
            $isoutput=0;
        }
    }
    else{
        if($isoutput==1){
            print outfile "$_\n";
        }
    }
}
close fasfile;
close outfile;
