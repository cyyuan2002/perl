#!/usr/bin/perl
use strict;

my $FasFile=shift;
open(my $fh_fasfile,$FasFile) || die "Can't open file $FasFile\n";
my $seqN;
my $seq;

while(<$fh_fasfile>){
    if(/^>(\S+)/){
        if ($seqN ne ""){
            open(my $fh_outfile,">$seqN.fa");
            print $fh_outfile ">$seqN\n$seq";
            close $fh_outfile;
        }
        $seqN=$1;
        $seq="";
    }
    else{
        $seq.=$_;
    }
}
open(my $fh_outfile,">$seqN.fa");
print $fh_outfile ">$seqN\n$seq";
close $fh_outfile;
close $fh_fasfile;
exit (0);
