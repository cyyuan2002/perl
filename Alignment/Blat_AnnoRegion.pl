#!/usr/bin/perl

#Bf_V2_1	14853	26377	-	6	570,307,188,573,140,82,	14853,16700,17617,19132,20272,26295,	jgi|Brafl1|85312|fgenesh2_pg.scaffold_140000064	gene_1
#Bf_V2_1	27745	33076	+	4	670,137,181,152,	27745,32290,32604,32924,	jgi|Brafl1|85313|fgenesh2_pg.scaffold_140000065	gene_2
#Bf_V2_1	37385	62028	-	6	152,171,240,236,190,250,	37385,45555,47340,53422,53885,61778,	jgi|Brafl1|85314|fgenesh2_pg.scaffold_140000066	gene_3
#Bf_V2_1	50392	51066	-	3	318,12,204,	50392,50726,50862,	jgi|Brafl1|208656|e_gw.36.51.1	gene_4
#Bf_V2_1	63661	64698	+	5	74,17,93,175,451,	63661,63735,63811,63909,64247,	jgi|Brafl1|72102|fgenesh2_pg.scaffold_36000157,jgi|Brafl1|85315|fgenesh2_pg.scaffold_140000067	gene_5
#Bf_V2_1	69205	76978	-	4	137,128,87,295,	69205,74724,75778,76683,	jgi|Brafl1|72101|fgenesh2_pg.scaffold_36000156,jgi|Brafl1|85316|fgenesh2_pg.scaffold_140000068	gene_6

use strict;

my ($parAnnofile,$parFlanklen,$parOutfile)=@ARGV;
if(@ARGV<2){
    print "Usage:$0 Anno_file Flank_length [Output_file]\n";
    exit;
}

if($parOutfile eq ""){
    $parOutfile=$parAnnofile.".rgn";
}

open (my $outfile,">$parOutfile");
open(my $annofile,$parAnnofile) || die "Can't open file: $parAnnofile\n";
while(<$annofile>){
    chomp();
    my @lines=split(/\t/,$_);
    my $left=$lines[2]-$parFlanklen;
    my $right=$lines[3]+$parFlanklen;
    if($left<1){
	$left=1;
    }
    if($right>$lines[1]){
	$right=$lines[1];
    }
    print $outfile "$lines[9]\t$lines[0]\t\+\t$left\t$right\n";
}
close $annofile;
close $outfile;
exit(0);
