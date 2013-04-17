#!/usr/bin/perl
use strict;
use File::Temp;
use File::Basename;

my ($BDFile,$DLFile)=@ARGV;
die "Usage:$0 <Breakdancer_CTX> <Delly_CTX> \n" if(@ARGV < 2);

die "Can't open file $BDFile\n" if(!-e $BDFile);
die "Can't open file $DLFile\n" if(!-e $DLFile);


my $BD_FILTER_SCORE=90;
my $DL_FILTER_SCORE=20;
my $DL_FILTER_COUNT=3;
my $MERGE_LENGTH=300;

my $Tempfile=File::Temp::tempnam(".","Temp");
$Tempfile=basename($Tempfile);
my $BDTemp1=$Tempfile.".BD1";
my $BDTemp2=$Tempfile.".BD2";

my $DLTemp1=$Tempfile.".DL1";
my $DLTemp2=$Tempfile.".DL2";

open(my $fh_bdtemp1,">$BDTemp1");
open(my $fh_bdfile, "$BDFile") || die "Can't open file $BDFile\n";
while(<$fh_bdfile>){
    next if(/^#/);
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[8] >= $BD_FILTER_SCORE){
        print $fh_bdtemp1 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\n";
    }
}
close $fh_bdfile;
close $fh_bdtemp1;

open(my $fh_dltemp1,">$DLTemp1");
open(my $fh_dlfile,"$DLFile");
while(<$fh_dlfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my @mapinfo=split(/\//,$lines[4]);
    my @scoreinfo=split(/\//,$lines[5]);
    my $scaf2=$mapinfo[0];
    my $site2=$mapinfo[1];
    my $mapscore=$scoreinfo[0];
    my $mapreads=$scoreinfo[1];
    if($mapscore >= $DL_FILTER_SCORE && $mapreads >= $DL_FILTER_COUNT){
        my @scaffold1=split(/\./,$lines[0]);
        my $chr1=$scaffold1[1];
        my @scaffold2=split(/\./,$scaf2);
        my $chr2=$scaffold2[1];
        my $site1=int(($lines[1]+$lines[2])/2);
        if($chr1 < $chr2){
            print $fh_dltemp1 "$lines[0]\t$site1\t$scaf2\t$site2\tCTX\n";
        }
        else{
            print $fh_dltemp1 "$scaf2\t$site2\t$lines[0]\t$site1\tCTX\n";
        }
    }
}
close $fh_dlfile;
close $fh_dltemp1;

my $sortBD=`cat $BDTemp1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BDTemp2`;
my $sortDL=`cat $DLTemp1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $DLTemp2`;

my %BD_SV;
my $lastchrom;
my $lastsite;
my $lastmatchchrom;
my $lastmatchsite;
open(my $fh_bdtemp2,"$BDTemp2");
while(<$fh_bdtemp2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[0] eq $lastchrom){
        if($lines[2] eq $lastmatchchrom){
            my $ss=$lastsite-$MERGE_LENGTH;
            my $se=$lastsite+$MERGE_LENGTH;
            my $es=$lastmatchsite-$MERGE_LENGTH;
            my $ee=$lastmatchsite+$MERGE_LENGTH;
            if($ss <= $lines[1] && $se >= $lines[1] && $es <= $lines[3] && $ee >=$lines[3]){
                next;
            }
        }
    }
    $lastchrom=$lines[0];
    $lastsite=$lines[1];
    $lastmatchchrom=$lines[2];
    $lastmatchsite=$lines[3];
    my %matchinfo;
    #print "$_\n";
    $matchinfo{'chr'}=$lines[2];
    $matchinfo{'site'}=$lines[3];
    $BD_SV{$lines[0]}->{$lines[1]}=\%matchinfo;
}
close $fh_bdtemp2;

open(my $fh_dltemp2,"$DLTemp2");
while(<$fh_dltemp2>){
    chomp();
    my @lines=split(/\t/,$_);
    my $chr=$lines[0];
    if(exists($BD_SV{$chr})){
        my %svsites=%{$BD_SV{$chr}};
        foreach my $site(sort {$a<=>$b} keys %svsites ){
            my $ss=$lines[1]-$MERGE_LENGTH;
            my $se=$lines[1]+$MERGE_LENGTH;
            if($ss <= $site && $se >= $site){
                my %svinfo=%{$svsites{$site}};
                if($lines[2] eq $svinfo{'chr'}){
                    my $es=$lines[3]-$MERGE_LENGTH;
                    my $ee=$lines[3]+$MERGE_LENGTH;
                    if($es <= $svinfo{'site'} && $ee >= $svinfo{'site'}){
                        print "$chr\t$site\t$svinfo{'chr'}\t$svinfo{'site'}\tCTX\n";
                    }
                }
            }
        }
    }
}

unlink($BDTemp1);
unlink($BDTemp2);
unlink($DLTemp1);
unlink($DLTemp2);

close $fh_dltemp2;
exit(1);
