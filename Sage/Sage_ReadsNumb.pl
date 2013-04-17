#!/usr/bin/perl
use strict;

my $file=shift;
open (file,$file) || die "Can't open file: $file\n";
my %sage;
my @reads=<file>;

for(my $i=0;$i<@reads;$i=$i+4){
    my $seq=$reads[$i+1];
    $seq=~s/\n//;
    if($sage{$seq} eq ""){
        $sage{$seq}=1;
    }
    else{
        $sage{$seq}=$sage{$seq}+1;
    }
}


foreach my $key(keys %sage){
    print "$key\t$sage{$key}\n";
}
