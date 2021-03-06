#!/usr/bin/perl
#===============================================================================
#
#         FILE: Seq_FqFilter.pl
#
#        USAGE: ./Seq_FqFilter.pl -a <fastq> -b <fastq> 
#
#  DESCRIPTION: Filter the low quality reads from illumia reads
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen 
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.1
#      CREATED: 01/15/2011 02:37:16 PM
#     REVISION: 05/10/2012
#===============================================================================

use strict;
use Compress::Zlib;
use Getopt::Long;

my $version="1.1 Yuan Chen 2012-05-10";
my %opts;
##Read parameters
#GetOptions(\%opts,"m=s","a=s","l:i","n:i","b:s","r:s","c:s","p:i","z:i","q:i","help");
GetOptions(\%opts,"m=s","a=s","l:i","n:i","b:s","c:s","p:i","z:i","q:i","help");

if((!defined $opts{m})||(!defined $opts{a})){
    &Usage();
}

my $parMode=$opts{m};
my $parFileA=$opts{a};
my $parFileB=$opts{b};
#my $parRemoveBs=(defined $opts{r}) ? $opts{r} : 40;
my $parRemoveNs=(defined $opts{n}) ? $opts{n} : 10;
my $parCutReads=(defined $opts{c}) ? $opts{c} : 0;
my $parRemovePolyA=(defined $opts{p}) ? $opts{p} : 0;
my $parOutZip=(defined $opts{z}) ? $opts{z} : 1;
my $parLowqFilter=(defined $opts{l}) ? $opts{l} : 1;
my $parQual=(defined $opts{q}) ? $opts{q} : 64;
my $parPolyALength=15;
my $parLowqFilterPercent=0.2;
my $parLowqFilterScore=20;


my $parFileMode; #input file mode
my @FileNA;
my @FileNB;
my $outFileA;
my $outFileB;
my $parCutStart;
my $parCutLength;
my $fhoutfileA;
my $fhoutfileB;
my $countBs=0;
my $countNs=0;
my $countPolyA=0;
my $countLowQ=0;
my $countAll=0;

my $rmCharacter=chr($parQual+2);

my %readsA;
my %readsB;



my $iserror=0;

if($parMode==2){
    if($parFileB eq ""){
	print stderr "Lack of fileb for pair-end mode\n";
	exit(1);
    }
}

##Check input file
&Filecheck($parFileA,0);
if($parMode==2){
    &Filecheck($parFileB,0);
}

if($iserror==1){
    exit(1);
}

##Autodetect input filetype, generate file name of output file;

@FileNA=split(/\./,$parFileA);
my $Filetype=$FileNA[scalar(@FileNA)-1];

if($Filetype eq "gz"){
    $parFileMode=1;
}
else{
    $parFileMode=0;
}


$outFileA=$FileNA[0];
for(my $i=1;$i<@FileNA-1;$i++){
    $outFileA="$outFileA.$FileNA[$i]";
}
$outFileA="$outFileA.cln";

if($parMode==2){
    @FileNB=split(/\./,$parFileB);
    $outFileB=$FileNB[0];
    for(my $i=1;$i<@FileNB-1;$i++){
        $outFileB="$outFileB.$FileNB[$i]";
    }
    $outFileB="$outFileB.cln";
}

if($parOutZip==1){#/*-/
    $outFileA="$outFileA.gz";
    $outFileB="$outFileB.gz";
    $fhoutfileA=gzopen($outFileA,"wb");
    $fhoutfileB=gzopen($outFileB,"wb");
}
else{
   open($fhoutfileA,">$outFileA");
   open($fhoutfileB,">$outFileB");
}

##build cut region
if($parCutReads!=0){
    my @Cut=split(/\:/,$parCutReads);
    $parCutStart=$Cut[0]-1;
    $parCutLength=$Cut[1];
}



##Reading Files;

if($parFileMode==1){ ## read gzip file
    if($parMode==1){ ##single-end;
	my $fhFileA=gzopen($parFileA,"rb");
        my $linecount;
        my %readsinfoA;
        while($fhFileA->gzreadline(my $lineA)){
            chomp($lineA);
            $linecount++;
            if($linecount==1){ $readsA{"N"}=$lineA;}
            elsif($linecount==2){ $readsA{"Q"}=$lineA;}
            elsif($linecount==3) {$readsA{"D"}=$lineA;}
            elsif($linecount==4){
                $readsA{"S"}=$lineA;
                my $fltreads=&readsflt();
		if($fltreads==1){
		    my $readsoutA=$readsA{"N"}."\n".$readsA{"Q"}."\n".$readsA{"D"}."\n".$readsA{"S"}."\n";
		    if($parOutZip==1){
			$fhoutfileA->gzwrite($readsoutA);
		    }
		    else{
			print $fhoutfileA "$readsoutA";
		    }
		}
		$linecount=0;
            }
        }
    }
    elsif($parMode==2){ ##pair-end;
	my $fhFileA=gzopen($parFileA,"rb");
        my $fhFileB=gzopen($parFileB,"rb");

	my $linecount;
        my %readsinfoA;
        while($fhFileA->gzreadline(my $lineA)){
	    my $lineB;
	    $fhFileB->gzreadline($lineB);
            chomp($lineA);
	    chomp($lineB);
            $linecount++;
            if($linecount==1){ $readsA{"N"}=$lineA;$readsB{"N"}=$lineB;}
            elsif($linecount==2){ $readsA{"Q"}=$lineA;$readsB{"Q"}=$lineB;}
            elsif($linecount==3) {$readsA{"D"}=$lineA;$readsB{"D"}=$lineB;}
            elsif($linecount==4){
                $readsA{"S"}=$lineA;
		$readsB{"S"}=$lineB;
		$countAll++;
                my $fltreads=&readsflt();

		if($fltreads==1){
		    my $readsoutA=$readsA{"N"}."\n".$readsA{"Q"}."\n".$readsA{"D"}."\n".$readsA{"S"}."\n";
		    my $readsoutB=$readsB{"N"}."\n".$readsB{"Q"}."\n".$readsB{"D"}."\n".$readsB{"S"}."\n";
		    if($parOutZip==1){
			$fhoutfileA->gzwrite($readsoutA);
			$fhoutfileB->gzwrite($readsoutB);
		    }
		    else{
			print $fhoutfileA "$readsoutA";
			print $fhoutfileB "$readsoutB";
		    }
		}
		$linecount=0;
            }
        }
	$fhFileA->gzclose();
	$fhFileB->gzclose();
	if($parOutZip==1){
	    $fhoutfileA->gzclose();
	    $fhoutfileB->gzclose();
	}
	else{
	    close $fhoutfileA;
	    close $fhoutfileB;
	}
	print "Total Reads: $countAll\nRemove reads contain $parRemoveNs Ns: $countNs\nRemove polyA reads: $countPolyA\nRemove low quality reads: $countLowQ\n";

	exit(0);

    }
    else{
	print stderr "Wrong parameters for mode, only 1 or 2 are allowed\n";
	exit(1);
    }
}
else {
    print stderr "Only support gzip file\n";
    exit(1);
}

sub readsflt{

    if($parMode==1){
        if($parCutReads==0){
	    my $stat=0;
#	    if($readsA{"Q"}=~/B{$parRemoveBs,}/){
#		$countBs++;
#		return 0;
#            }
	    if($readsA{"S"}=~/N{$parRemoveNs,}/){
		$countNs++;
		return 0;
            }
	    if($parRemovePolyA==1){
		my $headseq=substr($readsA{"Q"},0,$parPolyALength);
		if($headseq=~/A{$parPolyALength,}/ || $headseq=~/T{$parPolyALength,}/){
		    $countPolyA++;
		    return 0;
		}
	    }
	    if($parLowqFilter==1){
		my $lowscount=0;
		my @scoreA=split("",$readsA{"S"});
		for(my $i=0;$i<@scoreA;$i++){
		    my $score=ord($scoreA[$i])-$parQual;
		    if($score<$parLowqFilterScore){
			$lowscount++;
		    }
		}
		if($lowscount/length($readsA{"S"})>$parLowqFilterPercent){
		    $countLowQ++;
		    return 0;
		}
	    }
	    return 1;
	}
        else{
	    if($parRemovePolyA==1){
		my $headseq=substr($readsA{"Q"},0,$parPolyALength);
		if($headseq=~/A{$parPolyALength,}/ || $headseq=~/T{$parPolyALength,}/){
		    $countPolyA++;
		    return 0;
		}
	    }
            my $subAs=substr($readsA{"S"},$parCutStart,$parCutLength);
	    my $subAq=substr($readsA{"Q"},$parCutStart,$parCutLength);
            if($subAs=~/N{$parRemoveNs,}/){
		$countNs++;
                return 0;
            }
#	    if($subAq=~/B{$parRemoveBs,}/){
#		$countNs++;
#		return 0;
#            }
            else{
                $readsA{"S"}=$subAs;
		$readsA{"Q"}=$subAq;
		if($parLowqFilter==1){
		    my $lowscount=0;
		    my @scoreA=split("",$readsA{"S"});
		    for(my $i=0;$i<@scoreA;$i++){
			my $score=ord($scoreA[$i])-$parQual;
			if($score<$parLowqFilterScore){
			    $lowscount++;
			}
		    }
		    if($lowscount/length($readsA{"S"})>$parLowqFilterPercent){
			$countLowQ++;
			return 0;
		    }
		}
                return 1;
            }
        }
    }
    elsif($parMode==2){
        if($parCutReads==0){
#            if($readsA{"Q"}=~/B{$parRemoveBs,}/ || $readsB{"Q"}=~/B{$parRemoveBs,}/){
#		$countBs++;
#                return 0;
#            }
	    if($readsA{"S"}=~/N{$parRemoveNs,}/ || $readsB{"S"}=~/N{$parRemoveNs,}/){
		$countNs++;
                return 0;
            }
            if($parRemovePolyA==1){
		my $headseqA=substr($readsA{"Q"},0,$parPolyALength);
		my $offset=0-$parPolyALength;
		my $headseqB=substr($readsB{"Q"},$offset);
		if($headseqA=~/A{$parPolyALength,}/ || $headseqA=~/T{$parPolyALength,}/ || $headseqB=~/A{$parPolyALength,}/ || $headseqB=~/T{$parPolyALength,}/){
		    $countPolyA++;
		    return 0;
		}
	    }
	    if($parLowqFilter==1){
		my $lowscountA=0;
		my $lowscountB=0;
		my @scoreA=split("",$readsA{"S"});
		my @scoreB=split("",$readsB{"S"});
		for(my $i=0;$i<@scoreA;$i++){
		    my $scoreA=ord($scoreA[$i])-$parQual;
		    my $scoreB=ord($scoreB[$i])-$parQual;
		    if($scoreA<$parLowqFilterScore){
			$lowscountA++;
		    }
		    if($scoreB<$parLowqFilterScore){
			$lowscountB++;
		    }
		}
		if($lowscountA/length($readsA{"S"})>$parLowqFilterPercent || $lowscountB/length($readsA{"S"})>$parLowqFilterPercent){
		    $countLowQ++;
		    return 0;
		}
	    }

	    return 1;
        }
        else{
	    if($parRemovePolyA==1){
		my $headseqA=substr($readsA{"Q"},0,$parPolyALength);
		my $offset=0-$parPolyALength;
		my $headseqB=substr($readsB{"Q"},$offset);
		if($headseqA=~/A{$parPolyALength,}/ || $headseqA=~/T{$parPolyALength,}/ || $headseqB=~/A{$parPolyALength,}/ || $headseqB=~/T{$parPolyALength,}/){
		    $countPolyA++;
		    return 0;
		}
	    }
            my $subAs=substr($readsA{"S"},$parCutStart,$parCutLength);
	    my $subBs=substr($readsB{"S"},$parCutStart,$parCutLength);
	    my $subAq=substr($readsA{"Q"},$parCutStart,$parCutLength);
	    my $subBq=substr($readsB{"Q"},$parCutStart,$parCutLength);
#            if($subAq=~/B{$parRemoveBs,}/ || $subBq=~/B{$parRemoveBs,}/){
#		$countBs++;
#                return 0;
#            }
	    if($subAs=~/N{$parRemoveNs,}/ || $subBq=~/N{$parRemoveNs,}/){
		$countNs++;
		return 0;
	    }
            else{
                $readsA{"S"}=$subAs;
		$readsA{"Q"}=$subAq;
		$readsB{"S"}=$subBs;
		$readsB{"Q"}=$subBq;
		if($parLowqFilter==1){
		    my $lowscountA=0;
		    my $lowscountB=0;
		    my @scoreA=split("",$readsA{"S"});
		    my @scoreB=split("",$readsB{"S"});
		    for(my $i=0;$i<@scoreA;$i++){
			my $scoreA=ord($scoreA[$i])-$parQual;
			my $scoreB=ord($scoreB[$i])-$parQual;
			if($scoreA<$parLowqFilterScore){
			    $lowscountA++;
			}
			if($scoreB<$parLowqFilterScore){
			    $lowscountB++;
			}
		    }
		    if($lowscountA/length($readsA{"S"})>$parLowqFilterPercent || $lowscountB/length($readsA{"S"})>$parLowqFilterPercent){
			$countLowQ++;
			return 0;
		    }
		}
                return 1;
            }
        }
    }
}



sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "Directory $Filename already exists\n";
	    $iserror=1;
	}
    }
}



sub Usage(){
  print << "    Usage";

	This program is used to filter the reads in Fastq format. All the input files should be compressed by gzip. 

	Usage:  $0 (version $version)

	<options>
		-m     Program mode, 1 single-end File, 2 pair-end Files (must given)
		-a     File_name of Fastq file (must given)
		-b     File_name of Fastq file, only used when mode is 2 (must given when mode is 2)
		-n     Remove the reads with multiple Ns, (default 10)
		-c     Cut the reads length, start_site:length. (example: 2:75)
		-l     Remove the reads have more than 20% scores below 20 (0 off, 1 on), default 1
		-z     Output the clean reads in Zip format (0 off, 1 on), default 1
		-p     Remove PloyA in read (0 off, 1 on), default 0
		-q     Quality shift, default:64 <illumia 1.5>

    Usage

	exit(0);
};
