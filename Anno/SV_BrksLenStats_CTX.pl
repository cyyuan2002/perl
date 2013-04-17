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
my ($Filelst,$StrainOut,$OutFile)=@ARGV;
if(@ARGV<3){
    die "usage: <file_list> <Strain_Outfile> <Output_file>";
}
my $WINDOWSIZE=50;
my $TempFile1=$Filelst.".tmp1";
my $TempFile2=$Filelst.".tmp2";
my %allbreaks;

my $chromhead="supercont2.";
my $chromcount=14;

open(my $fh_filelst,"$Filelst");
open(my $fh_tmp1,">$TempFile1");
open(my $fh_strainout,">$StrainOut");

while(<$fh_filelst>){
    chomp();
    my $brkFile=$_;
    open(my $fh_brkfile,"$brkFile");
    my $lastchrom;
    my $svcount=0;
    my $svlength=0;
    my %chrsvcount;
    my %chrsvlength;
    while(<$fh_brkfile>){
        chomp();
        my @lines=split(/\t/,$_);
        if($lines[5] eq "pass"){
            print $fh_tmp1 "$_\n";
            if (exists($chrsvcount{$lines[0]})) {
                $chrsvcount{$lines[0]}++;
            }
            else{
                $chrsvcount{$lines[0]}=1;
            }
            if (exists($chrsvcount{$lines[2]})) {
                $chrsvcount{$lines[2]}++;
            }
            else{
                $chrsvcount{$lines[2]}=1;
            }
        }
    }
    close $fh_brkfile;

    $brkFile=~/aCneoH99\.r(\S+)\.\S+.brks/;
    my $strainName=$1;
    print $fh_strainout "$strainName";
    for(my $i=1;$i<=$chromcount;$i++){
        my $chromName=$chromhead.$i;
        if(exists($chrsvcount{$chromName})){
            print $fh_strainout "\t$chrsvcount{$chromName}";
        }
        else{
            print $fh_strainout "\t0";
        }
    }
    print $fh_strainout "\n";
}
close $fh_strainout;
close $fh_filelst;
close $fh_tmp1;

`cat $TempFile1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $TempFile2`;

open(my $fh_tmp2,"$TempFile2");
my @brks;
my $lastchrom;
while(<$fh_tmp2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[0] eq $lastchrom){
        my ($isoverlap,$index)=&checksite($lines[1],$lines[2],$lines[3]);
        if($isoverlap==0){
            my %info;
            $info{'posA'}=$lines[1];
            $info{'chr'}=$lines[2];
            $info{'posB'}=$lines[3];
            $info{'sv'}=$lines[4];
            $info{'count'}=1;
            push(@brks,\%info);
        }
        else{
            $brks[$index]->{'count'}++;
        }
    }
    else{
        my @tmparray=@brks;
        $allbreaks{$lastchrom}=\@tmparray;
        @brks=();
        $lastchrom=$lines[0];
        my %info;
        $info{'posA'}=$lines[1];
        $info{'chr'}=$lines[2];
        $info{'posB'}=$lines[3];
        $info{'sv'}=$lines[4];
        $info{'count'}=1;
        push(@brks,\%info);
    }
}
{
    my @tmparray=@brks;
    $allbreaks{$lastchrom}=\@tmparray;
}
close $fh_tmp2;

open(my $fh_rgn,">$OutFile");
foreach my $chrom (sort keys %allbreaks){
    my @brks=@{$allbreaks{$chrom}};
    my $SVcount;
    my $SVlength;
    for(my $i=0;$i<@brks;$i++){
        print $fh_rgn "$chrom\t",$brks[$i]->{'posA'},"\t",$brks[$i]->{'chr'},"\t",$brks[$i]->{'posB'},"\t",$brks[$i]->{'sv'},"\t",$brks[$i]->{'count'},"\n";
        $SVcount++;
    }
}
close $fh_rgn;

unlink ($TempFile1,$TempFile2);
exit(0);

sub checksite(){
    my ($posA,$chromB,$posB)=@_;
    for(my $i=0;$i<@brks;$i++){
        my %info=%{$brks[$i]};
        if ($chromB eq $info{'chr'}) {
            if($posA-$WINDOWSIZE <= $info{'posA'} && $posA+$WINDOWSIZE >=$info{'posA'}){
                if($posB-$WINDOWSIZE <= $info{'posB'} && $posB+$WINDOWSIZE >=$info{'posB'}){
                    return (1,$i);
                }
            }
        }
    }
    return 0;
}
