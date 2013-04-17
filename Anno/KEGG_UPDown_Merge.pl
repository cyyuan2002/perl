#!/usr/bin/perl
#===============================================================================
#
#         FILE: KEGG_UPDown_Merge.pl
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
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION:
#      CREATED:
#     REVISION:
#===============================================================================

use strict;
my ($upfile,$downfile,$KEGGfile)=@ARGV;
die "Usage: $0 <up_regulated_genes> <down_regulated_genes> <NCBIID_table>\n" if (@ARGV<3);

my %kegg;
open(my $fh_KEGGfile,$KEGGfile) || die "Can't open file $KEGGfile\n";
while(<$fh_KEGGfile>){
    chomp();
    my @lines=split(/\t/,$_);
    $kegg{$lines[0]}="$lines[1]\t$lines[2]";
}
close $fh_KEGGfile;

open(my $fh_upfile,$upfile) || die "Can't open file $upfile\n";
while(<$fh_upfile>){
    chomp();
    my @lines=split(/\t/,$_);
    print "$lines[0]\t$kegg{$lines[0]}\tU\n" if(exists($kegg{$lines[0]}));
}
close $fh_upfile;

open(my $fh_downfile,$downfile) || die "Can't open file $downfile\n";
while(<$fh_downfile>){
    chomp();
    my @lines=split(/\t/,$_);
    print "$lines[0]\t$kegg{$lines[0]}\tD\n" if(exists($kegg{$lines[0]}));
}
close $fh_downfile;

exit(1);
