#!/usr/bin/perl

#**************************************************************************
#This program is written by Chen Yuan (Email:chenyuan04@mails.gucas.ac.cn)
#Max-Planck Junior Scientist Group on Developmental Biology,Kunming Institute of Zoology,CAS
#Version 1.0
#File Created on 2005-11-15
#This program is used to draw maps based on blast result and output a text file for further analysis
#**************************************************************************

use GD; #For graphic output
use strict;
#-----------------Get user options-------------------------------------
use Getopt::Long;  

my %opts;
GetOptions(\%opts,"o:s","a:s","i:s","e:s","l:s","s:s","f:s","w:s","help");


if ((!defined $opts{i})|| (!defined $opts{o}) ) {
	Usage();
}

my $Input;
my $Output;
my $Expect;
my $Score;
my $SubLength;
my $Identity;
my $Filetype;
my $imWidth;

$Input = $opts{i};
$Output = $opts{o};
$Expect = (defined $opts{e}) ? $opts{e} : 5;
$Score = (defined $opts{s}) ? $opts{s} : 0;
$SubLength = (defined $opts{l}) ? $opts{l} : 30;
$Filetype = (defined $opts{f}) ? $opts{f} : "e";
$Identity = (defined $opts{a}) ? $opts{a} : 10;
$imWidth = (defined $opts{w}) ? $opts{w} : 1024;
#----------------------------------------------------------------------
#-----------------------Parameter List--------------------------------
print "Parameter List:\n";
print "Length\t-> $SubLength\n";
print "Expect\t-> $Expect\n";
print "Score\t-> $Score\n";
print "Identity\t-> $Identity\n";
print "File Type\t-> $Filetype\n";
print "Image Width\t-> $imWidth\n\n";
#---------------------------------------------------------------------

#---------------------------------variable list----------------------
open (FOpen,$Input) || die "cannot open file $Input\n";
print "Program is running...\n\n";

my %blastResult; #store all the result;
my $querySeq={}; #pointer of one query;
my @seqName; #array for subject query;
my $seqN;
my @seqAnnotation; 
my $seqAnno;
my $sLength;
my @seqLength; 
my $seqL;
my @seqScore;
my $seqS;
my @seqExpect;
my $expect;
my @identity1;
my $iden1;
my @identity2;
my $iden2;
my @idenpercent;
my $iden3;
my @qbegin;
my $qbeg;
my @qend;
my $qendl;
my @sbegin;
my $sbeg;
my @send;
my $sendl;
my @qBase;
my $qBa;
my @sBase;
my $sBa;

my $qStart;
my $sStart;
my $qSeq;
my $sSeq;
my $isSeq;
my $isQuery; 
my $isMatch;
my $seqCount;
my @blSeq;

$qStart=$sStart=$isSeq=$isQuery=$seqCount=$isMatch=0;
#-----------------------------------------------------------------------

#-----------------------read data for the file--------------------------
while(<FOpen>){
   if (/Query= \s?(\S+)/) {
	   #if there is new alignment sequence need to load, load it
   	   if ($isSeq==1) {
		   if($expect<=$Expect&&$seqS>=$Score&&$iden2>=$SubLength&&$iden3>=$Identity){
				push(@qBase,$qSeq);
				push(@sBase,$sSeq);
				push(@qbegin,$qbeg);
				push(@qend,$qendl);
				push(@sbegin,$sbeg);
				push(@send,$sendl);
				push(@seqScore,$seqS);	   
				push(@seqExpect,$expect);
	   			push(@identity1,$iden1);
				push(@identity2,$iden2);
				push(@idenpercent,$iden3);
				$seqCount++;				

				$isSeq=0;				
				$isMatch=1;
		   }
	   }
	   #if there is matched alignment in the subject sequence, load the suject sequence information
	   if ($isMatch==1) { 				
		   push(@blSeq,$seqCount);
	   	   push(@seqName,$seqN);
		   push(@seqLength,$seqL);
   		   push(@seqAnnotation,$seqAnno);
		   $isMatch=0
	   }
	   #if there is another query sequence, load the previous sequence information
	   if ($isQuery==1) { 
		   $querySeq->{seqName}=[@seqName];
		   $querySeq->{seqLength}=[@seqLength];
		   $querySeq->{seqScore}=[@seqScore];
		   $querySeq->{seqExpect}=[@seqExpect];
		   $querySeq->{identity1}=[@identity1];
		   $querySeq->{identity2}=[@identity2];
		   $querySeq->{idenpercent}=[@idenpercent];
		   $querySeq->{qBase}=[@qBase];
		   $querySeq->{sBase}=[@sBase];
		   $querySeq->{qbegin}=[@qbegin];
		   $querySeq->{qend}=[@qend];
		   $querySeq->{sbegin}=[@sbegin];
		   $querySeq->{send}=[@send];
		   $querySeq->{blSeq}=[@blSeq];
		   $querySeq->{seqAnnotation}=[@seqAnnotation];
		   $blastResult{$querySeq->{seqName}}=$querySeq;
		   $isQuery=0;
		   $querySeq={};
		   $seqCount=0;

		   @seqName=@seqLength=@seqScore=@seqExpect=@identity1=@identity2=@idenpercent=@qBase=@sBase=@qbegin=@qend=@sbegin=@send=@blSeq=@seqAnnotation=();
	   }
	   $querySeq->{queryName}=$1;
   }
   elsif(/\((\S+)\s+letters\)/){
	   my $q;
	   $q=$1;
	   $q=~s/,//;
	   $querySeq->{queryLength}=$q;
   }
   elsif(/^>(\S*)\s*(.*)/){ 
       if ($isSeq==1) {
   		   if($expect<=$Expect&&$seqS>=$Score&&$iden2>=$SubLength&&$iden3>=$Identity){

				push(@qBase,$qSeq);
				push(@sBase,$sSeq);
				push(@qbegin,$qbeg);
				push(@qend,$qendl);
				push(@sbegin,$sbeg);
				push(@send,$sendl);
				push(@seqScore,$seqS);	   
				push(@seqExpect,$expect);
	   			push(@identity1,$iden1);
				push(@identity2,$iden2);
				push(@idenpercent,$iden3);
				$seqCount++;

				$isSeq=0;
				$isMatch=1;
		   }
	   }	   
	   if ($isMatch==1) {				
		   push(@blSeq,$seqCount);	   
		   push(@seqName,$seqN);
		   push(@seqLength,$seqL);
		   push(@seqAnnotation,$seqAnno);
		   $isMatch=0
	   }

	   $seqN=$1;
	   $seqAnno=$2;
   }
   elsif(/Length\s?=\s?(\d+)/){
	   $seqL=$1;
   }
   elsif(/Score = (.+) bits.+Expect\S* =\s+(\S+)\s*/){
       if ($isSeq==1) {
   		   if($expect<=$Expect&&$seqS>=$Score&&$iden2>=$SubLength&&$iden3>=$Identity){

				push(@qBase,$qSeq);
				push(@sBase,$sSeq);
				push(@qbegin,$qbeg);
				push(@qend,$qendl);
				push(@sbegin,$sbeg);
				push(@send,$sendl);
				push(@seqScore,$seqS);	   
				push(@seqExpect,$expect);
	   			push(@identity1,$iden1);
				push(@identity2,$iden2);
				push(@idenpercent,$iden3);

				$seqCount++;
				$isSeq=0;
				$isMatch=1;
		   }
	   }
	   $seqS=$1;
	   $expect=$2;
	   $expect=~s/^e/1e/;
   }
   elsif (/Identities = (\d+)\/(\d+)\s+\((.{0,4})%\)/){

	   $iden1=$1;
	   $iden2=$2;
	   $iden3=$3;
	   $qStart=$sStart=0;
	   
   }
   elsif (/Query((\:\s+)|\s+)(\d+)\s*(\S+)\s+(\d+)/){

		if($qStart==0)
		{
			$qbeg=$3; 
			$qSeq=$4;
			
		}
		else { $qSeq .= $4;}
		$qendl=$5;
		$qStart=1;
	}
	elsif (/Sbjct((\:\s+)|\s+)(\d+)\s*(\S+)\s+(\d+)/){
		if($sStart==0)
		{
			$sbeg=$3;
			$sSeq=$4;
		}
		else {$sSeq.=$4;}
		$sendl=$5;
		$sStart=1;
		$isSeq=1;
		$isQuery=1;
	}	
}
if($expect<$Expect&&$seqS>$Score&&$iden2>$SubLength&&$iden3>=$Identity){

		   push(@qBase,$qSeq);
		   push(@sBase,$sSeq);
		   push(@qbegin,$qbeg);
		   push(@qend,$qendl);
		   push(@sbegin,$sbeg);
		   push(@send,$sendl);
		   push(@seqScore,$seqS);	   
		   push(@seqExpect,$expect);
	   	   push(@identity1,$iden1);
		   push(@identity2,$iden2);
		   push(@idenpercent,$iden3);
		   $seqCount++;

		   $isSeq=0;
		   $isMatch=1;
   }


if ($isMatch==1) {	
   push(@blSeq,$seqCount); 
   push(@seqName,$seqN);
   push(@seqLength,$seqL);
   push(@seqAnnotation,$seqAnno);
   $isMatch=0
}

$querySeq->{seqName}=[@seqName];
$querySeq->{seqLength}=[@seqLength];
$querySeq->{seqScore}=[@seqScore];
$querySeq->{seqExpect}=[@seqExpect];
$querySeq->{identity1}=[@identity1];
$querySeq->{identity2}=[@identity2];
$querySeq->{idenpercent}=[@idenpercent];
$querySeq->{qBase}=[@qBase];
$querySeq->{sBase}=[@sBase];
$querySeq->{qbegin}=[@qbegin];
$querySeq->{qend}=[@qend];
$querySeq->{sbegin}=[@sbegin];
$querySeq->{send}=[@send];
$querySeq->{blSeq}=[@blSeq];
$querySeq->{seqAnnotation}=[@seqAnnotation];
$blastResult{$querySeq->{seqName}}=$querySeq;

@seqName=@seqLength=@seqScore=@seqExpect=@identity1=@identity2=@idenpercent=@qBase=@sBase=@qbegin=@qend=@sbegin=@send=@blSeq=();
#-----------------------------------------------------------------------------------------------

my $key;
my $pointer;
my $qCount;
my $outFilename;

$qCount=0;
$outFilename=$Output;
my $outFile;

#call sub outtext to output a text file
&Outtext (\%blastResult,$Output,$Filetype);

#send each query sequence information to draw map
foreach $key (keys %blastResult) {
	$pointer=$blastResult{$key};
    $outFile=$outFilename.$qCount.".png";
	&draw($pointer,$outFile);
	$qCount++;
}
print "All jobs finished!\n";

#------------------------the sub used to draw picture-------------------
sub draw(){	
   my($point,$fileOut)=@_;
   my $keys;
   #----read data from the pointer-------------
   my @seqName;
   push(@seqName,$point->{queryName});
   foreach $keys (@{$point->{seqName}}) {
	   push(@seqName,$keys);
   }

   my @seq;
   push(@seq,$point->{queryLength});
   foreach $keys (@{$point->{seqLength}}) {
	   push (@seq,$keys);
   }

   my @sbegin;
   foreach $keys (@{$point->{sbegin}}) {
	   push (@sbegin,$keys);
   }

   my @send;
   foreach $keys (@{$point->{send}}) {
	   push (@send,$keys);
   }
   
   my @blSeq;
   foreach $keys (@{$point->{blSeq}}) {
	   push (@blSeq,$keys);
   }

   my @qbegin;
   foreach $keys (@{$point->{qbegin}}) {
	   push (@qbegin,$keys);
   }

   my @qend;
   foreach $keys (@{$point->{qend}}) {
	   push (@qend,$keys);
   }

   my $cne;
   my $subCount=1;
   my $blCount=0;
   my $i;
   for($i=0;$i<@sbegin;$i++){	   
	   if ($i==$blSeq[$blCount]) {
		   $subCount++;
		   $blCount++;
	   }
	   $cne->[$subCount][$i][0]=$sbegin[$i];
	   $cne->[$subCount][$i][1]=$send[$i];
   }

   for($i=0;$i<@qbegin;$i++){
	  $cne->[0][$i][0]=$qbegin[$i];
	  $cne->[0][$i][1]=$qend[$i];

   }

    my $seqnumber=@seqName;

    #create a canvas to draw
   	my $imHeight=150*$seqnumber;
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
    
    #find the longest sequence in the blast result
    my $maxlength=0;
    my $seqlength;
	foreach $seqlength (@seq) {
		if ($maxlength<$seqlength) {
			$maxlength=$seqlength;
		}
    }

	#caculate each sequence length to draw on the picture
	my $arrCount = 0;
    my @recWidth;
	while ($arrCount < @seq) {
		$recWidth[$arrCount]=$seq[$arrCount]/$maxlength*$recMaxwidth;
		
		$arrCount++;
    }
    
#-------------------draw sequence----------------------------------
    $arrCount=0;
	my $recHeight=10;
	my $recSpace=100;
	my $recLeftspace=5;
	my @recX1;
	my @recX2;
	my @recY1;
	my @recY2;
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




	while($arrCount < @recWidth){
		$recX1[$arrCount]=$recLeftspace;
		$recY1[$arrCount]=$arrCount*$recHeight+$recSpace*($arrCount+1);
		$recX2[$arrCount]=$recWidth[$arrCount]+$recLeftspace;
        $recY2[$arrCount]=$recHeight*($arrCount+1)+$recSpace*($arrCount+1);
        $seqTitle="$seqName[$arrCount] (Length = $seq[$arrCount] bp)";
		$im->string(gdSmallFont,$recX1[$arrCount],$recY1[$arrCount]-15,$seqTitle,$black);
		$im->rectangle($recX1[$arrCount],$recY1[$arrCount],$recX2[$arrCount],$recY2[$arrCount],$black);
		#print "rectangle $arrCount X1= ".$recX1[$arrCount]." Y1= ".$recY1[$arrCount]." X2= ".$recX2[$arrCount]." Y2= ".$recY2[$arrCount]."\n";
		$arrCount++;
	}

	#----------------------------------------------------------------------
	my $i;
	my $j;

    my @colors;
	my $g;
	my $b;
	my $r;
	#create random colors
	for ($i=0;$i<@qbegin;$i++) {
		$g=rand(255);
		$b=rand(255);
		$r=rand(255);
		$colors[$i]=$im->colorAllocate($r,$g,$b);
	}


	#----------------draw alignment blocks-------------------------------------------
	my $cneRecX1;
	my $cneRecX2;
	my $cneRecY1;
	my $cneRecY2;
	
    

	for ($j=0;$j<@qbegin;$j++) {
			$cneRecX1->[0][$j]=$cne->[0][$j][0]/$maxlength*$recMaxwidth+$recLeftspace;
			$cneRecY1->[0][$j]=$recY1[0];
            $cneRecX2->[0][$j]=$cne->[0][$j][1]/$maxlength*$recMaxwidth+$recLeftspace;
			$cneRecY2->[0][$j]=$recY2[0];
			$im->filledRectangle($cneRecX1->[0][$j],$cneRecY1->[0][$j],$cneRecX2->[0][$j],$cneRecY2->[0][$j],$colors[$j]);

	}
	
	$i=1;
	$blCount=0;
	for ($j=0;$j<@sbegin;$j++) {
		if($j==$blSeq[$blCount]){
			$i++;
			$blCount++;
		}
		$cneRecX1->[$i][$j]=$cne->[$i][$j][0]/$maxlength*$recMaxwidth+$recLeftspace;
		$cneRecY1->[$i][$j]=$recY1[$i];
        $cneRecX2->[$i][$j]=$cne->[$i][$j][1]/$maxlength*$recMaxwidth+$recLeftspace;
		$cneRecY2->[$i][$j]=$recY2[$i];
		if ($cneRecX1->[$i][$j]>$cneRecX2->[$i][$j]) {
			$im->filledRectangle($cneRecX2->[$i][$j],$cneRecY1->[$i][$j],$cneRecX1->[$i][$j],$cneRecY2->[$i][$j],$colors[$j]);
		}
		else{
			$im->filledRectangle($cneRecX1->[$i][$j],$cneRecY1->[$i][$j],$cneRecX2->[$i][$j],$cneRecY2->[$i][$j],$colors[$j]);
		}	
	}
	#---------------------------------------------------------------------------
    
	#-------------------------draw line-----------------------------------------------------------
	
	$i=1;
	$blCount=0;
	for ($j=0;$j<@sbegin;$j++) {
		if($j==$blSeq[$blCount]){
			$i++;
			$blCount++;
		}
		$im->line($cneRecX1->[0][$j],$cneRecY2->[0][$j],$cneRecX1->[$i][$j],$cneRecY1->[$i][$j],$colors[$j]);
	    $im->line($cneRecX2->[0][$j],$cneRecY2->[0][$j],$cneRecX2->[$i][$j],$cneRecY1->[$i][$j],$colors[$j]);
		$arrCount++;
	}


	open (OUT,">$fileOut");
	print "Output Image File: $fileOut\n\n";
	binmode OUT;
    print OUT $im->png;
	close OUT;
}

#----------the sub used to output a text file--------------------------
sub Outtext(){
	my ($sequence,$fileout,$filetype)=@_;
    my $point;
	my $blresult;
	my $filename;
	if ($filetype eq "e") {
		$filename=$fileout.".xls";
		open (FOut,">$filename");
		print FOut "Query-name\tLetter\tSbject-Name\tLength\tQueryS\tQueryE\tSbjctS\tSbjctE\tScore\tE-value\tOverlap/total\tIdentity\tQuery_sequence\tSbject_sequence\tSbjct_description\n"
	}
	elsif ($filetype eq "t") {
		$filename=$fileout.".txt";
		open (FOut,">$filename");
	}
    
	#read data from the pointer
	foreach $blresult (keys %{$sequence}) {
        $point=${$sequence}{$blresult};		
		my $keys;
		my $queryName;
		my $queryLength;

		$queryName=$point->{queryName};
		$queryLength=$point->{queryLength};

		my @seqName;
		foreach $keys (@{$point->{seqName}}) {
			push (@seqName,$keys);
		}

		my @seqAnnotation;
		foreach $keys (@{$point->{seqAnnotation}}) {
			push (@seqAnnotation,$keys);
		}

		my @seqLength;
		foreach $keys (@{$point->{seqLength}}) {
			push (@seqLength,$keys);
		}

		my @seqScore;
		foreach $keys (@{$point->{seqScore}}) {
			push (@seqScore,$keys);
		}

		my @seqExpect;
		foreach $keys (@{$point->{seqExpect}}) {
			push (@seqExpect,$keys);
		}
		
		my @identity1;
		foreach $keys (@{$point->{identity1}}) {
			push (@identity1,$keys);
		}

		my @identity2;
		foreach $keys (@{$point->{identity2}}) {
			push (@identity2,$keys);
		}

		my @idenpercent;
		foreach $keys (@{$point->{idenpercent}}) {
			push (@idenpercent,$keys);
		}

		my @qbegin;
		foreach $keys (@{$point->{qbegin}}) {
			push (@qbegin,$keys);
		}

		my @qend;
		foreach $keys (@{$point->{qend}}) {
			push (@qend,$keys);
		}

		my @sbegin;
		foreach $keys (@{$point->{sbegin}}) {
			push (@sbegin,$keys);
		}
		
		my @send;
		foreach $keys (@{$point->{send}}) {
			push (@send,$keys);
		}

		my @qbase;
		foreach $keys (@{$point->{qBase}}) {
			push (@qbase,$keys);
		}

		my @sbase;
		foreach $keys (@{$point->{sBase}}) {
			push (@sbase,$keys);
		}
        
		my @blSeq;
		foreach $keys (@{$point->{blSeq}}) {
			push (@blSeq,$keys);
		}
        #write data to output file
		my $i;
		my $seqCount=0;
		if ($filetype eq "e") {
			for ($i=0;$i<@qbegin;$i++) {
				if ($blSeq[$seqCount]==$i) {
					$seqCount++;
				}
				print FOut "$queryName\t$queryLength\t$seqName[$seqCount]\t$seqLength[$seqCount]\t$qbegin[$i]\t$qend[$i]\t$sbegin[$i]\t$send[$i]\t$seqScore[$i]\t$seqExpect[$i]\t$identity1[$i]/$identity2[$i]\t$idenpercent[$i]\t$qbase[$i]\t$sbase[$i]\t$seqAnnotation[$seqCount]\n";
			}
		}
		elsif ($filetype eq "t") {				
			print FOut "query : $queryName ($queryLength bp)\n\n";
			for ($i=0;$i<@qbegin;$i++) {
				if ($blSeq[$seqCount]==$i) {
					print FOut "\n>$seqName[$seqCount] ($seqLength[$seqCount]) $seqAnnotation[$seqCount]\n";
					$seqCount++;
				}
				print FOut ">query ($qbegin[$i]~$qend[$i])\n$qbase[$i]\n";
				print FOut ">sbjct ($sbegin[$i]~$send[$i])\n$sbase[$i]\n";
			}
			print FOut "\n";
		}
		
	}
	print "Output Text File: $filename\n\n";
	close FOut;

}

sub Usage() #help subprogram
{
    print << "    Usage";
    
	Version	  : 1.0
	Date	  : 2005-11-14
	Author    : Chen Yuan
	Work Group: Max-Planck Junior Scientist Group on Developmental Biology,
		        Kunming Institute of Zoology,CAS
	Contact   : Email:chenyuan04\@mails.gucas.ac.cn MSN:cy_yuan2002\@hotmail.com

    Function description :

		Analyse blast result, use pictures to show the relationship between the align sequences.

	Usage: $0 <options>

		-i     input of blast result file , must be given (String)
					
		-o     output of alignment sequence according with the condition , must be given (String)
		
		-l     the minimum match length , default value is 30
	
		-s     the minimum match score , default value is 0

		-e     the minimum expect value of the alignment , default value is 5

		-a     the minimun percentage of the alignment identity, default value is 10
			
		-f     the file type of the text file: e(for excel), t(for text) (default excel)

		-w     the width of output image file , defalut value is 1024

		-help  Show help

    Usage

	exit(0);
};		
