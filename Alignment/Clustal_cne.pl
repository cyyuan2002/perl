#!/usr/bin/perl

use Getopt::Long;

#version information------------------------------------------------
my $Version="1.1";
my $Date="2005-10-28";
my $Update="2005-10-29";
my $Author="Alvin Chen";
my $function="Find conserved element in genes in ClustalW result";
my $Contact="chenyuan04@mails.gucas.ac.cn";
my %align={};
#--------------------------------------------------------------------

#option get----------------------------------------------------------
my %opts;
GetOptions(\%opts,"i:s","o:s","w:s","s:s","help");

if ((!defined $opts{i})|| (!defined $opts{o})) {
	Usage();
}

$filename=$opts{i};
open (I,$filename)||die "can't open input file:$filename!";
$output=$opts{o};
open (O,">$output")||die "can't open output file:$output!";
my $linelength=0;
my $lastlength=60;
my $linenum=0;
my $window=$opts{w};
my $minScore=$opts{s};
if ($window==0) {
	$window=50;
}
if ($minScore==0) {
	$minScore=30;
}
print "input filename: ".$filename."\n";
print "output filename: ".$output."\n";
print "window size: ".$window."\n";
print "minimum match number: ".$minScore."\n\n"; 
#--------------------------------------------------------------------
#Read Clustal Data into hash-----------------------------------------
while(<I>){
	chomp;
	$linenum++;
	$line=$_;
	if(/(\S+)\s+(\S+)/){
		#if($1=~/[CLUSTAL|MUSCLE]/){
		if($1=~/\*/ || $1=~/CLUSTAL/ || $1=~/MUSCLE/){
			next;		
		}
		if(exists $align{$1}){
			$align{$1}.=$2;
			if (length($2)<$lastlength) {
				$lastlength=length($2);				
			}
    	}
		else{
			$align{$1}=$2;
			$linelength=length($2);
		}
	}	
}
close I;
#-----------------------------------------------------------------------

#Read compare stars from the file---------------------------------------
my $group = keys %align;
$group=$group+1;
my $linecount=0;
my $groupcount=1;
my $strStar="";
open (P,$filename);

while (<P>) {
   chomp;
   $linecount=$linecount+1;
   $number=$group*$groupcount+2;
   if ($linecount==$number) {
      $groupcount++;
	  my $string=$_;
	  if ($linecount==$linenum) {
		 $mlastlength=0-$lastlength;
		 $stars=substr($string,$mlastlength,$lastlength);
		 $strStar.=$stars;
	  }
	  else {
		 my $length=0-$linelength;
	     $stars=substr($string,$length,$linelength);
	     $strStar.=$stars; 	
      }
   }
}
close P;
#------------------------------------------------------------------------

#Count match score-------------------------------------------------------
my $numLength=length($strStar);
my $subStar="";
my $numScore=0;
my $increase=$window/2;
for ($i=0;$i<=$numLength-$window;$i+=$increase) {
    $subStar=substr($strStar,$i,$window);
	$subStar=~s/\s//g;	
	$curScore=length($subStar);
	if ($curScore>=$minScore) {
		$matchStart[$numScore]=$i;
		$matchScore[$numScore]=$curScore;
		$numScore++;
	}
}
#-------------------------------------------------------------------------

#output file---------------------------------------------------------------
my $matchnum=$numScore;
my $printline="";
print O "Window Size is: ".$window."\n";
print O "The minimum match number is: ".$minScore."\n\n";
my $nextloop=0;
my $matchnumber=0;
if ($matchnum==0) {
    print O "No sequence found under the given condition!\n";
}
else {
	for ($i=0;$i<$matchnum;$i++) {
		$matchS=$matchStart[$i];
		for($n=$i;$n<$matchnum;$n++){
	      $matchEnd=$matchStart[$n]+$window;
			$nextloop=$n;
			if ($matchEnd<$matchStart[$n+1]) {
				last;
			}		
		}	

		$matchRegion=$matchEnd-$matchS;
		$Stars=substr($strStar,$matchS,$matchRegion);
		$Stars=~s/\s//g;
		$matchnumber=length($Stars);
		$matchpercent=$matchnumber/$matchRegion;
		
		print O "Match region:".$matchStart[$i]."~".$matchEnd."\n";
        print O "Match number is: ".$matchnumber."\n";
		print O "Percentage of matches is: ".$matchnumber."/".$matchRegion." (".$matchpercent."%)\n\n";

		foreach (keys %align) {
			if ($align{$_} ne '') {
			   $printline=substr($align{$_},$matchS,$matchRegion);
			   #print O "matchS is: ".$matchS."\t matchRegion is: ".$matchRegion."\n";
			   #print O $printline."\n";
			   #print O $align{$_}."\n";
				print O "$_ \t $printline\n";
			}
		}
		$printline=substr($strStar,$matchS,$matchRegion);
		print O "marks \t\t $printline\n\n";
		$i=$nextloop;
	}
}
close O;
print "Job finished!";
#---------------------------------------------------------------------------

#help subprogram------------------------------------------------------------
sub Usage #help subprogram
{
    print << "    Usage";
    
	Version: $Version
	Date   : $Date
	Update : $Update
	Author : $Author
	Contact: $Contact
	
	Function description :
	
		$function

	Usage: $0 <options>

		-i     input of clustalw result file , must be given (String)
		 
			
		-o     output of alignment sequence according with the condition , must be given (String)

		
		-w     the window size of searching process , default value is 50

		
		-s     the minimum match score(must less than window size) , default value is 30


		-help  Show help

    Usage

	exit(0);
};		
#--------------------------------------------------------------------------
