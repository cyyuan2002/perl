#!/usr/bin/perl -w
use strict;

die "Usage:$0 <*_1.fq.reads.stat> <*_1.fq.clean.dup.stat> <JobName>\n" if @ARGV<3;

my $reads_stat =shift;
my $dup_stat =shift;
my $jobName =shift;

my @Stat;
open IN,$reads_stat or die "$!";
while(<IN>){
	next if !/^\d/;
	chomp;
	@Stat = split(/\t/);
}
close IN;

my $Duplicate;
my $Usable_reads = 0;
open IN,$dup_stat or die "$!";
while(<IN>){
	if (/Duplicate_reads:(\d+)/){
		$Duplicate = $1;
	}elsif(/Clean_reads:(\d+)/){
		$Usable_reads = $1*2;
	}else{
		next;
	}
}
close IN;

if ( $Usable_reads == 0 ){
	die "Format error: $dup_stat\n";
}

my $Read_len = $Stat[1];
my $Raw_reads = $Stat[0]*2 ;
my $Raw_bases = $Stat[2];
my $GC = $Stat[5]."_".$Stat[6];
my $Q20 = $Stat[3]."_".$Stat[4];
my $Ns_num = $Stat[7]/$Raw_reads * 2 * 100;
#print "$Ns_num\n";
my $Low_qual = $Stat[8]/$Raw_reads * 2 * 100;
my $Adapter = $Stat[9]/$Raw_reads * 2 * 100;
my $Small = $Stat[10]/$Raw_reads * 2 * 100;
my $Dup = $Duplicate/$Raw_reads * 2 * 100;
my $Usable_len = $Stat[11]."_".$Stat[12];
my $Usable_bases = $Usable_reads * ($Stat[11] + $Stat[12])/2;
$Raw_reads = $Raw_reads / 1000000;
$Raw_bases = $Raw_bases / 1000000;
$Usable_reads = $Usable_reads / 1000000;
$Usable_bases = $Usable_bases / 1000000;

open (my $out,">$jobName.stat");
print $out "Total Reads:\t$Raw_reads\n";
print $out "Read length:\t$Read_len\n";
print $out "GC content:\t$GC\n";
print $out "Q20 score:\t$Q20\n";
print $out "Ns Reads:\t$Ns_num\n";
print $out "Low quality Reads:\t$Low_qual\n";
print $out "Duplication:\t$Dup\n";
print $out "Adapter:\t$Adapter\n";
print $out "Small:\t$Small\n";
print $out "Usable lenth:\t$Usable_len\n";
print $out "Usable reads:\t$Usable_reads\n";
print $out "Usable bases:\t$Usable_bases\n";
close $out;

