#!/usr/bin/perl

##this program is used for find orthologous genes between two species by using the Bi-directory best hit.
##Input files are two blast results file, which is in m8 format.

use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"a:s","b:s","e:s","o:s","help");

if((!defined $opts{a})||(!defined $opts{b})){
	&Usage();
}

my $blsFileA=$opts{a};
my $blsFileB=$opts{b};
my $evalue=(defined $opts{e})?$opts{e}:1e-5;
my $outfile=(defined $opts{o})?$opts{o}:"$blsFileA.orth";

my %matchA;
my $lastID;

open (my $fh_blsFileA,$blsFileA) || die "Can't open file: $blsFileA\n";
while(<$fh_blsFileA>){
    my @lines=split(/\t/,$_);
    if($lastID ne $lines[0]){
	$lastID=$lines[0];
	if($lines[10]<=$evalue){
	    $matchA{$lines[0]}=$lines[1];
	}
    }
}
close $fh_blsFileA;

$lastID="";
open(my $fh_blsFileB,$blsFileB) || die "Can't open file: $blsFileB\n";
open(my $fh_outfile,">$outfile");
while(<$fh_blsFileB>){
    my @lines=split(/\t/,$_);
    if($lastID ne $lines[0]){
	$lastID=$lines[0];
	if($lines[0] eq $matchA{$lines[1]}){
	    print $fh_outfile "$lines[1]\t$lines[0]\n";
	}
    }
}
close $fh_blsFileB;
exit(1);

sub Usage(){
	print << "    Usage";

	Usage: $0 <options>

		-a     Blast fileA in m8 format

		-b     Blast fileB in m8 format

		-e     e-value cutoff (default: 1e-5)

		-o     Output file (default: fileA.orth)

		-help  Show help

    Usage

	exit(0);
}
