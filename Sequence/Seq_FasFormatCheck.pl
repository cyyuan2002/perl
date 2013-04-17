#!/usr/bin/perl
use strict;
use Getopt::Long;

my %opts;
my $version="1.0 Alvin Chen 2011-06-10";

GetOptions(\%opts,"i=s","help");
if((!defined $opts{i})){
    &Usage();
}

my $parInputfile=$opts{i};

open(my $fh_infile,$parInputfile) || die "Can't open file $parInputfile\n";
while(<$fh_infile>){
	chomp();
	if($_=~/>/g){
		if(pos($_)!=1){
			print "$_\n";
		}
	}
}
close $fh_infile;

sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-i     Input file for assembly

    Usage
	exit(1);
}