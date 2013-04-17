#!/usr/bin/perl 
use strict;

my $file=shift;
open(file,$file) || die "Can't open file $file\n";
while(<file>){
    my @lines=split(/\t/,$_);
    if($lines[0]/$lines[10] > 0.95){
        print $_;
    }
}
close file;
