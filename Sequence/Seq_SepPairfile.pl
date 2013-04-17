#!/usr/bin/perl
use strict;

my $parPairfile=shift;
my $fhPairfile;
my $fhOutfile1;
my $fhOutfile2;

open($fhPairfile,"$parPairfile") || die "Can't open file:$parPairfile\n";

my $outfile1="$parPairfile.1.fa";
my $outfile2="$parPairfile.2.fa";
open($fhOutfile1,">$outfile1");
open($fhOutfile2,">$outfile2");

my $linecount=0;
while(<$fhPairfile>){
	$linecount++;
	if($linecount<3){
		print $fhOutfile1 "$_";
	}
	elsif($linecount==3){
		print $fhOutfile2 "$_";
	}
	elsif($linecount==4){
		print $fhOutfile2 "$_";
		$linecount=0;
	}
}
close $fhOutfile1;
close $fhOutfile2;
close $fhPairfile;

exit(0);