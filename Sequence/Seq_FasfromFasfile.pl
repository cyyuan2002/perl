#!/usr/bin/perl

use strict;

my ($filename,$seqname,$seqout,$mode)=@ARGV;

if(@ARGV<2){
	print "Usage:$0 fastafile SeqNum OutFile Mode[0 (small file)|1 (big file)]\n";
	exit;
}

if($mode=""){
	$mode=0;
}

if($seqout eq ""){
	$seqout=$seqname.".fas";
}

open (seqname,$seqname) || die "Can't open file $seqname";
my @seqN;
while(<seqname>){
	chomp();
	my @seqname=split(/\t/,$_);
	my $seqN=$seqname[0];
	$seqN=~s/\s//g;
	push(@seqN,$seqname[0]);
}
close seqname;

open (file,$filename) || die "Can't open file $filename";
open (fileout,">$seqout");

if($mode==0){
	my $seqName;
	my $seq;
	my %seqs;
	while(<file>){
		if(/^>(\S+)/){
			if($seq ne ""){
				$seqs{$seqName}=$seq;
			}
			$seqName=$1;
			$seq="$_";
		}
		else{
			$seq=$seq.$_;
		}
	}
	$seqs{$seqName}=$seq;
	close file;
	for(my $i=0;$i<@seqN;$i++){
		if(exists($seqs{$seqN[$i]})){
			print fileout "$seqs{$seqN[$i]}";
		}
		else{
			print stderr "Can't find $seqN[$i]\n";
		}
	}

}
else{
	my $isfound=0;
	while(<file>){
	    chomp();
	    if(/^>(\S+)/){
		my $ismatch=0;
		for(my $i=0;$i<@seqN;$i++){
		    if($1 eq $seqN[$i]){
			$ismatch=1;
			print fileout "$_\n";
			last;
		    }
		}
		if($ismatch==1){
		    $isfound=1;
		}
		else{
		    $isfound=0;
		}
	    }
	    else{
		if($isfound==1){
		    print fileout "$_\n";
		}
	    }
	}
	close file;
}
exit(0);
