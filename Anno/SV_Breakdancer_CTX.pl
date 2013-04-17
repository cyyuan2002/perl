#!/usr/bin/perl
use strict;
use File::Basename;

my ($BDFile)=shift;
my $CTXFlanking=100;

my $minScore=90;
my $minReads=5;

my $Tempfile1=basename($BDFile).".tmp1";
my $Tempfile2=basename($BDFile).".tmp2";
my $Outfile=basename($BDFile).".CTX.mft";

open(my $fh_temp1,">$Tempfile1");
open(my $fh_bdfile,"$BDFile") || die "Can't open file: $BDFile\n";
while(<$fh_bdfile>){
    next if(/^#/);
    chomp();
    my @lines=split(/\t/,$_);
    next if($lines[6] ne "CTX");
    next if($lines[8] < $minScore || $lines[9] < $minReads);
    print $fh_temp1 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[8]\t$lines[9]\n";
}
close $fh_bdfile;
close $fh_temp1;

`cat $Tempfile1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $Tempfile2`;

my $lastchrom;
my $lastS;
my $lastTC;
my $lastTS;
my $lastReads;
my $lastScore;

open(my $fh_temp2,"$Tempfile2");
open(my $fh_out,">$Outfile");
while(<$fh_temp2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastchrom eq $lines[0] && $lastTC eq $lines[2]){
        if(($lines[1] >= $lastS-$CTXFlanking && $lines[1] <= $lastS+$CTXFlanking) && ($lines[3] >= $lastTS-$CTXFlanking && $lines[3] <= $lastTS+$CTXFlanking)){
            $lastScore=($lastReads*$lastScore+$lines[4]*$lines[5])/($lastReads+$lines[4]);
            $lastReads+=$lines[4];
            next;
        }
    }
    if($lastchrom ne ""){
        print $fh_out "$lastchrom\t$lastS\t$lastTC\t$lastTS\t$lastReads\t$lastScore\n";
    }
    $lastchrom = $lines[0];
    $lastS = $lines[1];
    $lastTC = $lines[2];
    $lastTS = $lines[3];
    $lastReads = $lines[4];
    $lastScore = $lines[5];
}
close $fh_temp2;
print $fh_out "$lastchrom\t$lastS\t$lastTC\t$lastTS\t$lastReads\t$lastScore\n";
close $fh_out;

unlink $Tempfile1;
unlink $Tempfile2;
