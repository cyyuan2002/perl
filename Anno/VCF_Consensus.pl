#!/usr/bin/perl
#===============================================================================
#
#         FILE: VCF_Consensus.pl
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
#      CREATED: Sep 4, 2012
#     REVISION:
#===============================================================================

use strict;
my ($vcf_mergefile) = shift;
die "Usage:$0 <Merged_VCF_file>\n" if($vcf_mergefile eq "");

open(my $fh_vcffile,$vcf_mergefile) || die "Can't open file $vcf_mergefile\n";
while(<$fh_vcffile>){
    chomp();
    if(/^##/) {print "$_\n";}
    else{
        my @lines=split(/\t/,$_);
        my $isskip=0;
        for(my $i=9;$i<@lines;$i++){
            if($lines[$i] eq "."){
                $isskip=1;
                last;
            }
        }
        print "$_\n" if($isskip==0);
    }
}
close $fh_vcffile;

exit(1);
