#!/usr/bin/perl
use strict;

my ($filein,$fileout)=@ARGV;
if(@ARGV<2){
    print stderr "Usage:$0 input_file output_file\n";
    exit(0);
}

my %dirs;
my %genes;
my ($fhin,$fhout);
open($fhin,$filein) || die "Can't open file $filein\n";
open($fhout,">$fileout");
while(<$fhin>){
    chomp();
    my @lines=split(/\t/);
    if(exists($genes{$lines[1]})){
	$genes{$lines[1]}.=",$lines[0]";
	$dirs{$lines[1]}.=",$lines[2]";
    }
    else{
	$genes{$lines[1]}=$lines[0];
	$dirs{$lines[1]}=$lines[2];
    }
}
close $fhin;

foreach my $key(keys %genes){
    print $fhout "$genes{$key}\t$key\t$dirs{$key}\n";
}
close $fhout;
exit(0);
