#!/usr/bin/perl

#**************************************************************************
#This program is written by Chen Yuan (Email:cheny04@post.kiz.ac.cn)
#Max-Planck Junior Scientist Group on Developmental Biology,Kunming Institute of Zoology,CAS
#Version 1.0
#File Created on 2006-04-16
#This program is used to draw listblast.pl
#**************************************************************************

use strict;
use GD;

my ($infile,$fileout)=@ARGV;

open(filein,$infile);
my @info=<filein>;
close filein;
my $i;
my (%queryname,%subname,$query,$sub,$querys,$querye,$subs,$sube);
my (%drawedsub,$drawedcount);
my $subcount;
my $maxlength;

for ($i=1;$i<@info;$i++) {
	my @lines=split(/\t/,$info[$i]);
	if ($subname{$lines[2]}==1) {
	}
	else{
		$subname{$lines[2]}=1;
		$subcount++;
		$maxlength=$lines[1];
		#print "$lines[8],$subcount\n";
	}
}
#-------------------------------------------------

my $imHeight=50*$subcount;
my $imWidth=1024;
my $im = new GD::Image($imWidth,$imHeight);
my $recMaxwidth=$imWidth*0.9;

#define colors used in program
my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);       
my $red = $im->colorAllocate(255,0,0);      
my $blue = $im->colorAllocate(0,0,255);
my $peachpuff = $im->colorAllocate(255,218,185);
my $blue2 = $im->colorAllocate(73,171,243);

$im->transparent($peachpuff);
$im->interlaced('true');


my $recHeight=10;
my $recSpace=20;
my $recLeftspace=5;
my $seqTitle;

#draw ruler; 
my $dLength;
$im->filledRectangle($recLeftspace,0,$recMaxwidth+$recLeftspace,5,$blue2);
$dLength=length($maxlength);
my $divideNum;
my $Num=substr($maxlength,0,1);

for (my $a=1;$a<$dLength;$a++) {
	$Num.=0;
}
my $dNum=$Num/10;
my $dX=$recLeftspace;

$im->line($dX,5,$dX,15,$black);
$im->string(gdSmallFont,$dX,15,"0",$black);

my $cNum;
$cNum=$maxlength/$dNum;
for ($a=1;$a<$cNum;$a++) {
	$dX=$dNum*$a/$maxlength*$recMaxwidth+$recLeftspace;
	$im->line($dX,5,$dX,15,$black);
	$im->string(gdSmallFont,$dX,15,$dNum*$a,$black);
}

if (($maxlength-$dNum*$cNum)>($dNum/2)) {	
	$dX=$recLeftspace+$recMaxwidth;
	$im->line($dX,5,$dX,15,$black);
	$im->string(gdSmallFont,$dX,15,$maxlength,$black);
}


#--------------------------------------------------


for ($i=1;$i<@info;$i++) {
	my @lines=split(/\t/,$info[$i]);
	$query=$lines[0];
	$querys=$lines[4];
	$querye=$lines[5];
	$sub=$lines[2];
	my $result=&draw($sub,$querys,$querye);
}

open (OUT,">$fileout");
print "Output Image File: $fileout\n\n";
binmode OUT;
print OUT $im->png;
close OUT;


sub draw{
	my($subname,$qs,$qe)=@_;
	if ($drawedsub{$subname}<1) {
		$drawedsub{$subname}=$drawedcount*$recHeight+$recSpace*($drawedcount+1)+40;
		$im->string(gdSmallFont,$recLeftspace,$drawedsub{$subname},$subname,$black);

		$drawedcount++;
	}
	my $width=$qe-$qs;
	my $recWidth=$width/$maxlength*$recMaxwidth;
	my ($recX1,$recX2,$recY1,$recY2);
	$recX1=$recLeftspace+$qs/$maxlength*$recMaxwidth;
	$recY1=$drawedsub{$subname};
	$recX2=$recX1+$recLeftspace+$recWidth;
    $recY2=$recY1+$recHeight;
	print "$subname,$qs,$qe\n";
	print "$recX1,$recY1,$recX2,$recY2\n";
	print "_________________________\n";
	$im->filledRectangle($recX1,$recY1,$recX2,$recY2,$red);
}

