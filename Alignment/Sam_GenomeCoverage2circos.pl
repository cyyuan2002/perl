#!/usr/bin/perl
#===============================================================================
#
#         FILE: Sam_GenomeCoverage.pl
#
#        USAGE:
#
#  DESCRIPTION:This program is used to caculate reads coverage on genome by the
#              results of samtools pile -f
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION:
#      CREATED:
#     REVISION:
#===============================================================================

use strict;
my ($infile,$windowsize,$outfile)=@ARGV;
if(@ARGV<1){
    print "Usage: $0 input_file window_size(default 1000) output_file\n";
    exit(0);
}

$windowsize||=1000;
$outfile||="$infile.covr";

my $coversummer;
my @avercovers;
my $kcounter;
my $lastcontig;
my %contigcovers;
open(my $fh_infile,"$infile") || die "$!\n";
while(<$fh_infile>){
    my @lines=split(/\t/,$_);
    if($lastcontig ne $lines[0]){
        if(@avercovers>0){
            my @avers=@avercovers;
            $contigcovers{$lastcontig}=\@avers;
        }
        @avercovers=();
        $lastcontig=$lines[0];
        $kcounter=0;
    }
    if(($lines[1]-$kcounter*$windowsize)>$windowsize){
        $kcounter++;
        my $avercov=$coversummer/$windowsize;
        push(@avercovers,$avercov);
        $coversummer=0;
    }
    $coversummer+=$lines[3];
}
close $fh_infile;
my $avercov=$coversummer/$windowsize;
push(@avercovers,$avercov);
$contigcovers{$lastcontig}=\@avercovers;

open (my $fh_out,">$outfile");
foreach my $chrom (sort keys %contigcovers){
    my @covers=@{$contigcovers{$chrom}};
    my $kcount=1;
    for(my $i=0;$i<@covers;$i++){
        my $end=$kcount+$windowsize-1;
        print $fh_out "$chrom $kcount $end $covers[$i]\n";
        $kcount+=$windowsize;
    }
}
close $fh_out;
exit(1);
