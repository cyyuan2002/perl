#!/usr/bin/perl
#===============================================================================
#
#         FILE:
#
#        USAGE:
#
#  DESCRIPTION: This program is used to deal with histo data
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
my $file=shift;
my $outfile="$file.cov";

open(my $fh_file,$file) || die "$!\n";
open(my $fh_out,">$outfile");
my $lastchrom;
my @pos;
my @scores;
while(<$fh_file>){
    chomp();
    my @lines=split(" ",$_);
    if($lastchrom ne $lines[0]){
        if($lastchrom ne ""){
           my @sortscore=sort {$a<=>$b} @scores;
           my $midpos=int((scalar(@sortscore))/2);
           print "$sortscore[$midpos]\n";
           for(my $i=0;$i<@scores;$i++){
            $scores[$i]=$scores[$i]-$sortscore[$midpos];
            print $fh_out "$pos[$i] $scores[$i]\n";
           }
        }
        @scores=();
        @pos=();
        $lastchrom=$lines[0];
    }
    my $score=pop(@lines);
    push(@scores,$score);
    push(@pos,join(" ",@lines));
}
my @sortscore=sort {$a<=>$b} @scores;
my $midpos=int((scalar(@sortscore))/2);
print "$sortscore[$midpos]\n";

for(my $i=0;$i<@scores;$i++){
    $scores[$i]=$scores[$i]-$sortscore[$midpos];
    print $fh_out "$pos[$i] $scores[$i]\n";
}
close $fh_file;
close $fh_out;
exit(1);
