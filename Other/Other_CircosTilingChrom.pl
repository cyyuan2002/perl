#!/usr/bin/perl
#===============================================================================
#
#         FILE: Other_CircosTilingChrom.pl
#
#        USAGE: This program is used to create chrom information for circos by the result of Mummer_Tiling.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.0
#      CREATED: 30/05/12
#     REVISION:
#===============================================================================

use strict;
my $infile=shift;
my $gapcolor="grey";
my $mapcolor="white";

open(my $fh_infile,$infile) || die "Can't open file $infile\n";
my @chromIDs;
my @chromlength;
my @Bands;
my %chromBands;
my $lastchrom;
my $gapstart=0;
my $lastchromlength;
my $count;
while(<$fh_infile>){
    chomp();
    if(/^>(\S+)\s\S+\s\S+\s(\S+)/){
        if($lastchrom ne ""){
            my $gapband="band $lastchrom gap$count gap$count $gapstart $lastchromlength $gapcolor";
            push(@Bands,$gapband);
            my @tempBands=@Bands;
            $chromBands{$lastchrom}=\@tempBands;
        }
        $lastchrom=$1;
        push(@chromIDs,$1);
        push(@chromlength,$2);
        @Bands=();
        $gapstart=0;
        $lastchromlength=$2;
        $count=0;
    }
    else{
        my @lines=split(/\t/,$_);
        my $gapend=$lines[2]-1;
        if($gapend>$gapstart){
            my $gapband="band $lastchrom gap$count gap$count $gapstart $gapend $gapcolor";
            push(@Bands,$gapband);
        }
        my $mapband="band $lastchrom map$count map$count $lines[2] $lines[3] $mapcolor";
        push(@Bands,$mapband);
        $gapstart=$lines[3]+1;
        $count++;
    }
}
my @tempBands=@Bands;
$chromBands{$lastchrom}=\@tempBands;
close $fh_infile;

my $iscolor=0;
$iscolor=1 if(@chromIDs<24);
for(my $i=0;$i<@chromIDs;$i++){
    print "chr - $chromIDs[$i] $chromIDs[$i] 0 $chromlength[$i]";
    print " chr",$i+1 if($iscolor);
    print "\n";
}

foreach my $chrom(sort keys %chromBands){
    my @Bandsinfo=@{$chromBands{$chrom}};
    foreach my $band(@Bandsinfo){
        print "$band\n";
    }
}

exit(1);
