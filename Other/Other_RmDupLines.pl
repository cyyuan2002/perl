#!/usr/bin/perl
##This script is used to remove the lines with duplicates, it will leave the first line with the same id

use strict;
my ($infile,$colnum)=@ARGV;
die "Usage:$0 <input file> <column number>\n" if(@ARGV<2);
my %ids;
my $rmcount=0;

open(my $fh_infile,$infile) || die "Can't open file $infile\n";
while(<$fh_infile>){
    chomp();
    my @lines=split(/\t/);
    if(!exists($ids{$lines[$colnum-1]})){
        print "$_\n";
        $ids{$lines[$colnum-1]}=1;
    }
    else{
        $rmcount++;
    }
}
print stderr "$rmcount lines removed from the data\n";
exit(1);

