#!/usr/bin/perl

use strict;

my ($filename,$seqname,$mode)=@ARGV;

if(@ARGV<2){
	print "Usage:$0 fastafile seqname mode(0|1)\n";
	exit;
}

if($mode eq ""){
	$mode=0;
}

open (file,$filename) || die "Can't open file $filename";
my $isfound=0;
while(<file>){
	chomp();
	if(/^>/){
		if($isfound==1 && $mode==0){
			exit;
		}
		if(/$seqname/){
			$isfound=1;
			print "$_\n";
		}
		else{
			$isfound=0;
		}
	}
	else{
		if($isfound==1){
			print "$_\n";
		}
	}
}

