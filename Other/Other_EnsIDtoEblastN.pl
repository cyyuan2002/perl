#!/usr/bin/perl
use strict;

my($parEnsidfile,$parEblastNfile)=@ARGV;
open(my $ensfile,$parEnsidfile) || die "Can't open file $parEnsidfile\n";
my %geneName;
my %geneDes;

while(<$ensfile>){
	chomp();
	my @lines=split(/\t/,$_);
	$geneName{$lines[0]}=$lines[1];
	$geneDes{$lines[0]}=$lines[2];
}
close $ensfile;

open(my $blastfile,$parEblastNfile) || die "Can't open file $parEblastNfile\n";
while(<$blastfile>){
	chomp();
	my @lines=split(/\t/,$_);
	print "$lines[0]\t$geneName{$lines[11]}\t$geneDes{$lines[11]}\t$lines[11]\n";
}
close $blastfile;

exit(0);