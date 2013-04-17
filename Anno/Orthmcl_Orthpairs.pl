#!/usr/bin/perl
#===============================================================================
#
#         FILE: Orthmcl_Orthpairs.pl
#
#        USAGE: 
#
#  DESCRIPTION:This program is used for searching orthologous pairs by given
#               species names from the results of orthmcl
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

my ($orthfile,$speA,$speB)=@ARGV;
die "Usage: $0 orthmcl_group.txt speciesA speciesB" if(@ARGV<3);
open(my $fh_orthfile,$orthfile);
while(<$fh_orthfile>){
    chomp();
    if(/$speA/ && /$speB/){
        my @lines=split(" ",$_);
        my @speAs;
        my @speBs;
        foreach my $gid(@lines){
            my @gname=split(/\|/,$gid);
            push(@speAs,$gname[1]) if($gname[0] eq $speA);
            push(@speBs,$gname[1]) if($gname[0] eq $speB);
        }
        print join(",",@speAs),"\t",join(",",@speBs),"\n";
    }
}
close $fh_orthfile;
exit(1);