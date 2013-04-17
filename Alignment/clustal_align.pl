#!/usr/bin/perl

use Getopt::Long;

my $Version="1.0";
my $Date="2005-10-29";
my $Author="Alvin Chen";
my $function="Export ClustalW result in line";
my $Contact="chenyuan04@mails.gucas.ac.cn";
my %align={};

my %opts;
GetOptions(\%opts,"i:s","o:s","help");
if ((!defined $opts{i})|| (!defined $opts{o})) {
	Usage();
}

$filename=$opts{i};
open (I,$filename)||die "can't open input file:$filename!";
$output=$opts{o};
open (O,">$output")||die "can't open output file:$output!";
my $lastlength=60;
my $linenum=0;

print "input filename: ".$filename."\n";
print "output filename: ".$output."\n\n";

while(<I>){
	chomp;
	$linenum++;
	$line=$_;
	if(/(\S+)\s+(\S+)/){
		if($1=~/CLUSTAL|\*/){
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
		}
	}	
}
close I;

my $group = keys %align;
$group=$group+1;
my $strlength=60;
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
		 #$align{"star"}.=$stars;
	  }
	  else {
		 my $length=0-$strlength;
	     $stars=substr($string,$length,$strlength);
	     $strStar.=$stars; 	
      }
   }
}
close P;

foreach  (keys %align) {
	if ($align{$_} ne '') {
		print O "$_ \t $align{$_}\n";
	}
}

close O;
print "File $output created!";
sub Usage #help subprogram
{
    print << "    Usage";
    
	Version: $Version
	Date   : $Date
	Author : $Author
	Contact: $Contact
	
	Function description :
	
		$function

	Usage: $0 <options>

		-i     input of clustalw result file , must be given (String)
		 
			
		-o     output of alignment sequence filter out gaps , must be given (String)

		
		-help  Show help

    Usage

	exit(0);
};		
