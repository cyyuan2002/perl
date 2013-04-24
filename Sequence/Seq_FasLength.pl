#!/usr/bin/perl
use strict;
my $filein=shift;

open(filein,$filein);
my $seqN;
my $seqLen=0;
while(<filein>){
    chomp();
    if(/^>(\S+)/){
        if ($seqLen > 0) {
            print "$seqN\t$seqLen\n";
        }
        
        $seqN=$1;
        $seqLen=0;
    }
    else{
        $seqLen+=length($_);
    }
}
print "$seqN\t$seqLen\n";
close filein;

exit(0);