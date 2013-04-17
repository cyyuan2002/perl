#!/usr/bin/perl
use strict;
my %transName;
my %transDes;
my %geneName;
my %geneDes;

my ($parTransAnno,$parGeneAnno,$parGeneID)=@ARGV;
open(my $transanno,$parTransAnno) || die "Can't open file: $parTransAnno\n";
while(<$transanno>){
	chomp();
	my @lines=split(/\t/,$_);
	$transName{$lines[0]}=$lines[1];
	$transDes{$lines[0]}=$lines[2];
}
close $transanno;

open(my $geneanno,$parGeneAnno) || die "Can't open file: $parGeneAnno\n";
while(<$geneanno>){
	chomp();
	my @lines=split(/\t/,$_);
	$geneName{$lines[0]}=$lines[1];
	$geneDes{$lines[0]}=$lines[2];
}
close $geneanno;

open (my $geneid,$parGeneID) || die "Can't open file: $parGeneID\n";
<$geneid>;
while(<$geneid>){
	chomp();
	my $id=$_;
	$id=~s/\"//g;
	if(exists ($geneName{$id})){
		print "$id\t$geneName{$id}\t$geneDes{$id}\n";
	}
	elsif(exists ($transName{$id})){
		print "$id\t$transName{$id}\t$transDes{$id}\n";
	}
	else{
		print "$id\tNoDes\n";
	}
}
close $geneid;

exit(0);
