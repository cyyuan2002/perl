#!/usr/bin/perl
use strict;
my ($filename,$chrominfo)=@ARGV;

my %chroms;
my $windowsize=5000;

open(my $fh_chrominfo,$chrominfo) || die "Can't open file $chrominfo";
my @chromName;
while(<$fh_chrominfo>){
    chomp();
    my @lines=split(/\t/,$_);
    $chroms{$lines[0]}=$lines[1];
    push(@chromName,$lines[0]);
}
close $fh_chrominfo;

my %Snps;

my $lastchrom;
my %chromsnps;
open(my $fh_vcffile,$filename);
while(<$fh_vcffile>){
    chomp();
    next if(/^#/);
    my @lines=split(/\t/,$_);
    my $chrom=$lines[0];
    if($lastchrom ne $chrom){
        if($lastchrom ne ""){
            $Snps{$lastchrom}=\%chromsnps;
        }
        $lastchrom=$lines[0];
        %chromsnps={};
    }

    my $site=int($lines[1]/$windowsize);
    if(!exists($chromsnps{$site})){
        $chromsnps{$site}=1;
    }
    else{
        $chromsnps{$site}++;
    }
}
$Snps{$lastchrom}=\%chromsnps;
close $fh_vcffile;

foreach my $chrom(@chromName){
    if(exists($Snps{$chrom})){
        my %rep=%{$Snps{$chrom}};
        my $loop=int($chroms{$chrom}/$windowsize)+1;
        for(my $i=0;$i<$loop;$i++){
            if(exists($rep{$i})){
                my $start=$i*$windowsize;
                my $end=($i+1)*$windowsize-1;
                print "$chrom $start $end $rep{$i}\n";
            }
            else{
                my $start=$i*$windowsize;
                my $end=($i+1)*$windowsize-1;
                print "$chrom $start $end 0\n";
            }
        }
    }
}
