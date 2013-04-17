#!/usr/bin/perl
use strict;

my ($filein,$fileout,$mode)=@ARGV;
if(@ARGV<3){
    print stderr "Usage:$0 input_file output_file mode(0|1)\n";
    exit(0);
}

my ($fhin,$fhout);
open($fhin,$filein) || die "Can't open file $filein\n";
open($fhout,">$fileout");

my $QueryN;
my $QueryS;
while(<$fhin>){
    my @lines=split(/\t/,$_);
    if($QueryN ne $lines[0]){
	print $fhout "$_";
	$QueryN=$lines[0];
	$QueryS=$lines[7];
    }
    else{
	if($mode==1) {
	    if($QueryS<=$lines[7]){
	        print $fhout "$_";
	    }
	}
    }
}
close $fhin;
close $fhout;
