#!/usr/bin/perl
use strict;

my ($Filelst,$OutFile,$SoftNum)=@ARGV;
die "Usage:$0 <File_list> <Output_file> <Software_Num>" if(@ARGV<2);

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

open(my $fh_file,"$Filelst");
open(my $fh_out,">$OutFile");
while(<$fh_file>){
    chomp();
    my $brkFile=$_;
    my %softcount;
    my %softpasscount;
    my $totalcount=0;
    my $passcount=0;
    my $totallength=0;
    my $SVtype;
    open(my $fh_brkfile,$brkFile) || die "Can't find file: $brkFile\n";
    while(<$fh_brkfile>){
        chomp();
        my @lines=split(/\t/,$_);
        my @soft=uniq(split(/:/,$lines[6]));
        $SVtype=$lines[4];
        my $softnum=scalar(@soft);
        $totalcount++;
        if($lines[5] eq "pass"){
            $passcount++;
            if(!exists($softpasscount{$softnum})){
                $softpasscount{$softnum}=1;
            }
            else{
                $softpasscount{$softnum}++;
            }
            if($lines[4] ne "CTX"){
                $totallength+=abs($lines[3]-$lines[1]+1);
            }
        }
        if(!exists($softpasscount{$softnum})){
            $softcount{$softnum}=1;
        }
        else{
            $softcount{$softnum}++;
        }
    }
    $brkFile=~/aCneoH99\.r(\S+)\.\S+.brks/;
    my $strainName=$1;
    
    if($SVtype ne "CTX"){
        print $fh_out "$strainName\t$totalcount\t$passcount\t$totallength\t";
    }
    else{
        print $fh_out "$strainName\t$totalcount\t$passcount\t";
    }
    for(my $i=$SoftNum;$i>0;$i--){
        if(!exists($softpasscount{$i})){
            print $fh_out "0/0(0%)\t";
        }
        else{
            my $percent=$softpasscount{$i}/$softcount{$i}*100;
            my $per=sprintf("%.2f",$percent);
            print $fh_out "$softpasscount{$i}/$softcount{$i}($per%)\t"; 
        }
    }
    print $fh_out "\n";
}
close $fh_out;
exit(0);