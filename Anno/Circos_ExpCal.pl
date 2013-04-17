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

#-----Format of the input file
##chromA pos1 pos2 exp_value
#supercont2.1 19966 20490 1
#supercont2.1 46817 309447 1
#supercont2.1 102790 103523 1
#supercont2.1 108768 108831 2
#supercont2.1 113467 113658 4
#------------------------

#----Format of the chrom info
#chromA length
#supercont2.1	2291499
#supercont2.2	1621675
#supercont2.3	1575141
#supercont2.4	1084805
#supercont2.5	1814975
#--------------

use strict;

my ($expFile,$chrominfo,$windowsize)=@ARGV;
die "Usage:$0 <Exp_File> <Chrom_Info> [Windowsize, default: 1000]\n" if (@ARGV<3);

my %chroms;
$windowsize ||=1000;

my @chromName;
open(my $fh_chrominfo,$chrominfo) || die "Can't open file $chrominfo";
while(<$fh_chrominfo>){
    chomp();
    my @lines=split(/\t/,$_);
    $chroms{$lines[0]}=$lines[1];
    push(@chromName,$lines[0]);
}
close $fh_chrominfo;

my %ChromRep;

open(my $fh_expfile,$expFile);
while(<$fh_expfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my $chrom=$lines[0];
    my $start=int($lines[1]/$windowsize);
    my $end=int($lines[3]/$windowsize)+1;
    for(my $i=$start;$i<$end;$i++){
        if(!exists($ChromRep{$chrom})){
            my %rep;
            $rep{$i}=$lines[5];
            $ChromRep{$chrom}=\%rep;
        }
        else{
            my %rep=%{$ChromRep{$chrom}};
            if(!exists($rep{$i})){
                $rep{$i}=$lines[5];
            }
            else{
                $rep{$i}+=$lines[5];
            }
            $ChromRep{$chrom}=\%rep;
        }
    }
}
close $fh_expfile;

foreach my $chrom(@chromName){
    if(exists($ChromRep{$chrom})){
        my %rep=%{$ChromRep{$chrom}};
        my $loop=int($chroms{$chrom}/$windowsize)+1;
        for(my $i=0;$i<$loop;$i++){
            my $start=$i*$windowsize;
            my $end=($i+1)*$windowsize-1;
            if(exists($rep{$i})){
                print "$chrom $start $end $rep{$i}\n";
            }
            else{
                print "$chrom $start $end 0\n";
            }
        }
    }
}

exit(0);

