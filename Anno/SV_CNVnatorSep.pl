#!/bin/perl
##This script is used to seperate DEL/DUP from CNVnator vcf results 

use strict;

my $file=shift;
my $DUP_out="$file.dup";
my $DEL_out="$file.del";

open(my $fh_outdup,">$DUP_out");
open(my $fh_outdel,">$DEL_out");
open(my $fh_filein,"$file");
while(<$fh_filein>){
	if(/^#/){
		print $fh_outdel $_;
		print $fh_outdup $_;
	}
	else{
		if(/\<DEL\>/){
			print $fh_outdel $_;
		}
		elsif(/\<DUP\>/){
			print $fh_outdup $_;
		}
		else{
			print $_;
		}
	}
}
close $fh_filein;
close $fh_outdup;
close $fh_outdel;

exit(0);