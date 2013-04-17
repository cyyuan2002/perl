#!/usr/bin/perl
use strict;

my $file=shift;
my $count=0;

open(my $fh_file,"$file");
while(<$fh_file>){
    chomp();
    my @lines=split(/\t/,$_);
    my $end1=$lines[1]+1;
    my $end2=$lines[3]+1;
    print "tranloc$count $lines[0] $lines[1] $end1\n";
    print "tranloc$count $lines[2] $lines[3] $end2\n";
    $count++;
}
close $fh_file;

exit(0);