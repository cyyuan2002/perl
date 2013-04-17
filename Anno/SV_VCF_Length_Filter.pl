#!/usr/bin/perl
use strict;

my ($VCF_File,$Min_Length,$SV_Type)=@ARGV;
die "Useage: <VCF_File> <Min_Length> <SV_Type>" if(@ARGV < 3);

my $Outfile="$VCF_File.flt";
open(my $fh_vcf,"$VCF_File") || die "Can't find file: $VCF_File\n";
open(my $fh_out,">$Outfile");
while(<$fh_vcf>){
    chomp();
    if(/^#/){
        print $fh_out "$_\n";
        next;
    }
    chomp();
    my @lines=split(/\t/,$_);
    my @infos=split(/;/,$lines[7]);
    my $endsite=0;
    for my $info(@infos){
        if($info=~/END=(\d+)/){
            $endsite=$1;
            last;
        }
    }
    if($endsite-$lines[1] >= $Min_Length){
        if(/SVTYPE=$SV_Type/){
            print $fh_out "$_\n";
        }
    }
}
close $fh_vcf;

