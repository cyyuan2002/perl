#!/usr/bin/perl
use strict;

my $file=shift;
open(my $fh_file,$file) || die "Can't open file\n";
my $lastid;
my $ids;
my %gos;
while(<$fh_file>){
	chomp();
	my @lines=split(/\t/,$_);
	my $gid=$lines[0];
	$gid=~/(\w+_\d+)T\d/;
	$gid=$1;
	if($lastid ne $gid){
		if($lastid ne ""){
			$gos{$lastid}=$ids;
		}
		$lastid=$gid;
		$ids=$lines[4];
	}
	else{
		$ids.=", $lines[4]";
	}
}
$gos{$lastid}=$ids;

foreach my $key(keys %gos){
	print "$key\t$gos{$key}\n";
}