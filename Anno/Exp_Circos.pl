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

my ($ExpFile,$ChromFile,$GTFFile)=@ARGV;
die "Usage: <Express_File> <Chrom_Length> <GTF_File>" if(@ARGV<3);

my $windowssize=50000;

open(my $fh_ExpFile,"$ExpFile") || die "Can't open file: $ExpFile\n";
my %gene_exp;
while(<$fh_ExpFile>){
    chomp();
    $_=~/\"(\S+)\"\t(\S+)/;
    $gene_exp{$1}=$2;
}
close $fh_ExpFile;

my %chroms;
open(my $fh_ChromFile,"$ChromFile") || die "Can't open file: $ChromFile\n";
while(<$fh_ChromFile>){
    chomp();
    my @lines=split(/\t/,$_);
    $chroms{$lines[0]}=$lines[1];
}
close $fh_ChromFile;

my %Geneannos;
open(my $fh_GTFFile,"$GTFFile") || die "Can't open file: $GTFFile\n";
my $geneS=0;
my $geneE=0;
my $geneName="";
my $Chrom="";
my @genes;
while(<$fh_GTFFile>){
    chomp();
    if(/start_codon/){
        my @lines=split(/\t/,$_);
        if($Chrom ne $lines[0]){
            if($Chrom ne ""){
                my @temparray;
                if(exists($Geneannos{$Chrom})){
                    @temparray=@{$Geneannos{$Chrom}};
                    push(@temparray,@genes);
                }
                else{
                    @temparray=@genes;
                }
                $Geneannos{$Chrom}=\@temparray;
            }
            $Chrom=$lines[0];
            @genes=();
        }
        $geneS=$lines[3];
        if($lines[8]=~/gene_id \"(\S+)\"/){
            $geneName=$1;
        }
    }
    elsif(/stop_codon/){
        my @lines=split(/\t/,$_);
        $geneE=$lines[4];
        my %geneinfo;
        $geneinfo{'id'}=$geneName;
        if($geneS<$geneE){
            $geneinfo{'s'}=$geneS;
            $geneinfo{'e'}=$geneE;
        }
        else{
            $geneinfo{'s'}=$geneE;
            $geneinfo{'e'}=$geneS;
        }
        push(@genes,\%geneinfo);
    }
}
{
    my @temparray=@genes;
    $Geneannos{$Chrom}=\@temparray;
    @genes=();
}
close $fh_GTFFile;

foreach my $chrom (sort keys %chroms){
    my $chromlength=$chroms{$chrom};
    my $posS=1;
    my $posE=$posS+$windowssize;
    while($posS < $chromlength){
        if($posE > $chromlength){
            $posE=$chromlength;
        }
        my @genes=@{$Geneannos{$chrom}};
        my $genecount=0;
        my $totalexp=0;
        foreach my $refgeneinfo (@genes){
            my %geneinfo=%{$refgeneinfo};
            next if($geneinfo{'e'} < $posS || $geneinfo{'s'} > $posE);
            if(exists($gene_exp{$geneinfo{'id'}})){
                $genecount++;
                $totalexp+=$gene_exp{$geneinfo{'id'}};
            }
        }
        if($genecount > 0){
            my $aveexp=$totalexp/$genecount;
            print "$chrom $posS $posE $aveexp\n";
        }
        else{
            print "$chrom $posS $posE 0\n";
        }
        $posS+=$windowssize;
        $posE+=$windowssize;
    }
}

exit(0);
