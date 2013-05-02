#!/usr/bin/perl
use strict;

my %GOs;
my $GOTable=shift;
open(my $fh_in,$GOTable) or die "Can't open file: $GOTable";
while (<$fh_in>) {
    chomp();
    my @lines=split(/\t/,$_);
    next if (/\tEC:(\S+)/);
    if (!exists($GOs{$lines[0]})) {
        $GOs{$lines[0]}=$lines[1];
    }
    else{
        $GOs{$lines[0]}.=", $lines[1]";
    }
}
close $fh_in;

foreach my $key(sort keys %GOs){
    print "$key\t$GOs{$key}\n";
}
