#!/usr/bin/perl

#===============================================================================
#
#         FILE: SV_DellyMerge.pl
#
#        USAGE:
#
#  DESCRIPTION: This script is used to merge and filter delly output
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

my ($Filein,$Mode,$SupportingReads)=@ARGV;
if(@ARGV < 2){
    die "Usage $0 <Delly_File> <Mode: 0 for DEL & DUP, 1 for INV,2 for CTX> [Supporting_Reads (default:5)]\n";
}
my $Fileout="$Filein.mft";

my $OutFlanking=50;
my $InFlanking=100;
my $CTXFlanking=100;
my $QualityFilter=20;

my $minlength=50;
my $mergeoverlap=0.5;
#my $maxlength=1000000;

$SupportingReads ||=5;

my $Tempfile1="$Filein.tmp";
my $Tempfile2="$Filein.tmp2";
open(my $fh_filein,"$Filein") || die "Can't open file: $Filein\n";
open (my $fh_outtemp1,">$Tempfile1");

while(<$fh_filein>){
    if(/\>\S+\</){
        print $fh_outtemp1 "$_";
    }
}
close $fh_filein;
close $fh_outtemp1;

if($Mode==2){
    `cat $Tempfile1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $Tempfile2`;
}
else{
    `cat $Tempfile1 | sort -k 1,1 -k 2,2n -k 3,3n > $Tempfile2`;
}

my $lastchrom;
my $lastS;
my $lastE;
my $lastReads;
my $lastScore;
my $lastTC;
my $lastTS;
my $lastCover;

open (my $fh_temp2, "$Tempfile2");
open (my $fh_out, ">$Fileout");
while (<$fh_temp2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($Mode==0){ ##DEL,DUP
        next if($lines[4]<$SupportingReads);
        next if($lines[5] <$QualityFilter);
        next if($lines[2]-$lines[1] < $minlength);
        if($lastchrom eq $lines[0]){
            if($lines[1] >= $lastE){
                my $length=$lastE-$lastS;
                print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastReads\t$lastScore\n";
                $lastchrom=$lines[0];
                $lastS=$lines[1];
                $lastE=$lines[2];
                $lastReads=$lines[4];
                $lastScore=$lines[5];
            }
            if($lines[3] <= $lastE){
                $lastScore=($lastReads*$lastScore+$lines[4]*$lines[5])/($lastReads+$lines[4]);
                $lastReads+=$lines[4];
            }
            else{
                my $overlap=$lastE-$lines[1];
                my $lastlength=$lastE-$lastS;
                my $newlength=$lines[3]-$lines[1];
                if($overlap/$lastlength >= $mergeoverlap && $overlap/$newlength >= $mergeoverlap){
                    $lastE=$lines[2];
                    $lastReads+=$lines[4];
                }
                else{
                    my $length=$lastE-$lastS;
                    print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastReads\t$lastScore\n";
                    $lastchrom=$lines[0];
                    $lastS=$lines[1];
                    $lastE=$lines[3];
                    $lastReads=$lines[4];
                    $lastScore=$lines[5];
                }
            }
        }
        if($lastchrom ne $lines[0]){
            if($lastchrom ne ""){
                my $length=$lastE-$lastS;
                print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastReads\t$lastScore\n";
            }
            $lastchrom = $lines[0];
            $lastS = $lines[1];
            $lastE = $lines[2];
            $lastReads = $lines[4];
            $lastScore = $lines[5];
        }
    }
    elsif($Mode==1){ #INV
        next if($lines[5] <$QualityFilter);
        if($lastchrom eq $lines[0]){
            if(($lines[1] >= $lastS-$InFlanking && $lines[1] <= $lastS+$OutFlanking) && ($lines[2] >= $lastE-$InFlanking && $lines[2] <= $lastE+$OutFlanking) ){
                if($lines[2] > $lastE){
                    $lastE=$lines[2];
                }
                $lastScore=($lastReads*$lastScore+$lines[4]*$lines[5])/($lastReads+$lines[4]);
                $lastCover=($lastCover+$lines[4])/2;
                next;
            }
        }
        if($lastchrom ne ""){
            if($lastScore >= $QualityFilter){
                my $length=$lastE-$lastS;
                print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastCover\t$lastScore\n";
            }
        }
        $lastchrom = $lines[0];
        $lastS = $lines[1];
        $lastE = $lines[2];
        $lastCover = $lines[4];
        $lastScore = $lines[5];
    }
    else{ ##CTX
        next if($lines[5] < $QualityFilter || $lines[4] < $SupportingReads);
        if($lastchrom eq $lines[0] && $lastTC eq $lines[2]){
            if(($lines[1] >= $lastS-$CTXFlanking && $lines[1] <= $lastS+$CTXFlanking) && ($lines[3] >= $lastTS-$CTXFlanking && $lines[3] <= $lastTS+$CTXFlanking)){
                $lastScore=($lastReads*$lastScore+$lines[4]*$lines[5])/($lastReads+$lines[4]);
                $lastReads+=$lines[4];
                next;
            }
        }
        if($lastchrom ne ""){
            if($lastReads >= $SupportingReads && $lastScore >= $QualityFilter){
                print $fh_out "$lastchrom\t$lastS\t$lastTC\t$lastTS\t$lastReads\t$lastScore\n";
            }
        }
        $lastchrom = $lines[0];
        $lastS = $lines[1];
        $lastTC = $lines[2];
        $lastTS = $lines[3];
        $lastReads = $lines[4];
        $lastScore = $lines[5];
    }
}
if($Mode==0){
    if($lastReads >= $SupportingReads && $lastScore >= $QualityFilter){
        my $length=$lastE-$lastS;
        print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastReads\t$lastScore\n";
    }
}
elsif($Mode==1){
    if($lastScore >= $QualityFilter){
        my $length=$lastE-$lastS;
        print $fh_out "$lastchrom\t$lastS\t$lastE\t$length\t$lastCover\t$lastScore\n";
    }
}
else{
    if($lastReads >= $SupportingReads && $lastScore >= $QualityFilter){
        print $fh_out "$lastchrom\t$lastS\t$lastTC\t$lastTS\t$lastReads\t$lastScore\n";
    }
}

close $fh_temp2;
close $fh_out;

unlink $Tempfile1;
unlink $Tempfile2;

exit(0);
