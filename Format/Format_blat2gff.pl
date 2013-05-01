#!/usr/bin/perl

#===============================================================================
#
#         FILE: 
#
#        USAGE:
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 
#     REVISION:
#===============================================================================

use strict;

my ($pslFile,$gffFile) =@ARGV;
if (@ARGV < 2) {
    die "Usage: <pslFile> <outGFF>\n";
}

open(my $fh_out,">"."$gffFile");
print $fh_out "##gff-version 3\n";

open(my $fh_psl,"$pslFile") or die "Can't open file:$pslFile\n";
while (<$fh_psl>) {
    next unless(/^\d/);
    chomp();
    my @lines= split(/\t/,$_);
    my $geneName=$lines[9];
    my $dir=$lines[8];
    my $chrom=$lines[13];
    my $chrom_s=$lines[15];
    my $chrom_e=$lines[16];
    my @blocklength=split(/,/,$lines[18]);
    my @blockstarts=split(/,/,$lines[20]);
    print $fh_out "$chrom\tBLAT\tgene\t$chrom_s\t$chrom_e\t.\t$dir\t.\tID=$geneName;Name=$geneName\n";
    print $fh_out "$chrom\tBLAT\tmRNA\t$chrom_s\t$chrom_e\t.\t$dir\t.\tID=$geneName","T0;Parent=$geneName\n";
    
    if ($dir eq "-") {
        @blocklength=reverse(@blocklength);
        @blockstarts=reverse(@blockstarts);
    }
    for(my $i=0;$i<@blockstarts;$i++){
        my $end=$blockstarts[$i]+$blocklength[$i];
        print $fh_out "$chrom\tBLAT\texon\t$blockstarts[$i]\t$end\t.\t$dir\t.\tID=exon:$geneName","T0:$i;Parent=$geneName","T0\n";
    }
}
close $fh_psl;
close $fh_out;
exit(0);
