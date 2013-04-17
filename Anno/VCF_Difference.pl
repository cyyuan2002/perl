#!/usr/bin/perl
#===============================================================================
#
#         FILE: VCF_Difference.pl
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
#      CREATED: Sep 5, 2012
#     REVISION:
#===============================================================================


use strict;
my ($vcf_mergefile, $groupAnum, $groupBnum) = @ARGV;
die "Usage:$0 <Merged_VCF_file> <GroupA_num> <GroupB_num>\n" if(@ARGV<3);

open(my $fh_vcffile,$vcf_mergefile) || die "Can't open file $vcf_mergefile\n";
my $lastAcolumn=9+$groupAnum-1;
my @groupAcolumn=(9..$lastAcolumn);
my $firstBcolumn=$lastAcolumn+1;
my $lastBcolumn=$firstBcolumn+$groupBnum-1;
my @groupBcolumn=($firstBcolumn..$lastBcolumn);
while(<$fh_vcffile>){
    chomp();
    if(/^##/) {print "$_\n";}
    else{
        my @lines=split(/\t/,$_);
        my $isskip=1;
        

    }
}
close $fh_vcffile;

exit(1);
