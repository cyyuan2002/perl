#!/usr/bin/perl
#===============================================================================
#
#         FILE: Mummer_CalCoords.pl
#
#        USAGE:Mummer_CalCoords.pl -c <coords_file> -r <reference_seqs> -q <query_seqs> -f [1:List, 2:Tabular]
#
#  DESCRIPTION: This script is used to calculate the identity and coverage of
#               the coord file, which is obtained by MUMmer software package
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.1
#      CREATED: 05-20-2012
#     REVISION: 10-1-2012
#===============================================================================

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"c=s","r=s","q=s","f:i","help");
if(!defined($opts{c}) || !defined($opts{r}) || !defined($opts{q})){
    die "Usage: $0 -c <coords_file> -r <reference_seqs> -q <query_seqs> -f [1:List, 2:Tabular]\n";
}

my $coordfile=$opts{c};
my $reffile=$opts{r};
my $queryfile=$opts{q};
my $outstyle=$opts{f};

$outstyle ||=1;

my %reflength;
my $reftotallength;
my $seqID;
my $seq;
my $seqlength;
my @refseqID;

open(my $fh_reffile,$reffile) || die "$!\n";
while(<$fh_reffile>){
    chomp();
    if(/^>(\S+)/){
	if($seq ne ""){
	    $seqlength=length($seq);
	    $reflength{$seqID}=$seqlength;
	    $reftotallength+=$seqlength;
	}
	$seqID=$1;
	$seq="";
	push(@refseqID,$seqID);
    }
    else{
	$seq.=$_;
    }
}
close $fh_reffile;
$seqlength=length($seq);
$reflength{$seqID}=$seqlength;
$reftotallength+=$seqlength;
$seqID="";
$seq="";

my $querytotallength;
my %querylength;

open(my $fh_queryfile,$queryfile) || die "$!\n";
while(<$fh_queryfile>){
    chomp();
    if(/^>(\S+)/){
	if($seq ne ""){
	    $seqlength=length($seq);
	    $querylength{$seqID}=$seqlength;
	    $querytotallength+=$seqlength;
	}
	$seqID=$1;
	$seq="";
    }
    else{
	$seq.=$_;
    }
}
close $fh_queryfile;
$seqlength=length($seq);
$querylength{$seqID}=$seqlength;
$querytotallength+=$seqlength;

my %refcover;
my %refidentity;
my %refmapID;
my %querymapID;

my $scafID;
my $scaf_e;
my $over_len;
my $scaf_len;
my @identity;
my @mappinglength;
my $totalcoverlength;
my $totalsimlength;
my $totalmappinglength;
my %mapids;

open(my $fh_coordfile,$coordfile) || die "$!\n";  ##coverage 要计算overlap，identity使用mapping的长度Xidentity累加除以mapping总长
while(<$fh_coordfile>){
    chomp();
    my @lines=split(/\t/,$_);
    if($scafID ne $lines[11]){ ##new scaf
	if($scafID ne ""){ #skip first record
	    #my $coverage=$over_len/$scaf_len*100;
	    $totalcoverlength+=$over_len;
	    $refcover{$scafID}=$over_len;
	    my $maplength;
	    my $similarlength;
	    for(my $i=0;$i<@identity;$i++){
		my $sim=$mappinglength[$i]*$identity[$i]/100;
		$similarlength+=$sim;
		$maplength+=$mappinglength[$i];
	    }
	    my $simall=$similarlength/$maplength*100;
	    $refidentity{$scafID}=$simall;
	    $totalsimlength+=$similarlength;
	    $totalmappinglength+=$maplength;
	    $refmapID{$scafID}=join(",",(keys %mapids));
	}
	$scafID=$lines[11];
	$scaf_len=$lines[7];
	$over_len=$lines[1]-$lines[0]+1;
	$scaf_e=$lines[1];
	@identity=();
	@mappinglength=();
	push(@identity,$lines[6]);
	push(@mappinglength,$over_len);
	%mapids=();
	$mapids{$lines[12]}=1;
	$querymapID{$lines[12]}=1 if(!exists($querymapID{$lines[12]}));
    }
    else{
	if($lines[0]<=$scaf_e){
	    if($lines[1]>=$scaf_e){
		$over_len=$over_len+($lines[1]-$scaf_e);
		$scaf_e=$lines[1];
		my $maplen=$lines[1]-$lines[0]+1;
		push(@identity,$lines[6]);
		push(@mappinglength,$maplen);
		$mapids{$lines[12]}=1 if(!exists($mapids{$lines[12]}));
		$querymapID{$lines[12]}=1 if(!exists($querymapID{$lines[12]}));
	    }
	}
	else{
	    $over_len=$over_len+($lines[1]-$lines[0]);
	    $scaf_e=$lines[1];
	    my $maplen=$lines[1]-$lines[0]+1;
	    push(@identity,$lines[6]);
	    push(@mappinglength,$maplen);
	    $mapids{$lines[12]}=1 if(!exists($mapids{$lines[12]}));
	    $querymapID{$lines[12]}=1 if(!exists($querymapID{$lines[12]}));
	}
    }
}

#my $coverage=$over_len/$scaf_len*100;
$refcover{$scafID}=$over_len;
my $maplength;
my $similarlength;
for(my $i=0;$i<@identity;$i++){
    my $sim=$mappinglength[$i]*$identity[$i]/100;
    $similarlength+=$sim;
    $maplength+=$mappinglength[$i];
}
my $simall=$similarlength/$maplength*100;
$refidentity{$scafID}=$simall;
$totalcoverlength+=$over_len;
$totalsimlength+=$similarlength;
$totalmappinglength+=$maplength;
close $fh_coordfile;

my $totalcoverper=sprintf("%.2f",$totalcoverlength/$reftotallength*100);
my $totalqueryper=sprintf("%.2f",$totalcoverlength/$querytotallength*100);
my $totalsimper=sprintf("%.2f",$totalsimlength/$totalmappinglength*100);

my @unmappedref;
foreach my $key (sort {$a<=>$b} keys %reflength){
    if(!exists($refcover{$key})){
	push (@unmappedref,$key);
    }
}
my $unmaplengthref;
foreach my $key(@unmappedref){
    $unmaplengthref+=$reflength{$key};
}
my $unmapreflengthper=sprintf("%.2f",$unmaplengthref/$reftotallength*100);

my @unmappedqr;
foreach my $key (sort {$a<=>$b} keys %querylength){
    if(!exists($querymapID{$key})){
	push (@unmappedqr,$key);
    }
}
my $unmapqrlength;
foreach my $key(@unmappedqr){
    $unmapqrlength+=$querylength{$key};
}
my $unmapqrlengthper=sprintf("%.2f",$unmapqrlength/$querytotallength*100);

my $mappedqrnum=scalar(keys %querylength)-scalar(@unmappedqr);

if($outstyle==1){
    print "Referenece: $reffile\n";
    print "Query: $queryfile\n";
    print "Number of sequences: ",scalar(keys %reflength)," ",scalar(keys %querylength),"\n";
    print "Total length: $reftotallength $querytotallength\n";
    print "Total alignment length: $totalcoverlength\n";
    print "Total alignment percentage: $totalcoverper\% $totalqueryper\%\n";
    print "Similarity: $totalsimper\%\n";
    print "Mapped sequences: ",scalar(keys %refcover),"/",scalar(keys %reflength)," ",$mappedqrnum,"/",scalar(keys %querylength),"\n";
    
    
    print "\nUnmapped reference sequence number: ",scalar(@unmappedref),"\n";
    print "Unmapped reference length: $unmaplengthref $unmapreflengthper\%\n";
    print "Sequence IDs: ";
    scalar(@unmappedref) < 1 ? print "\n": print join(",",@unmappedref),"\n";
    
    print "\nUnmapped query sequence number: ",scalar(@unmappedqr),"\n";
    print "Unmapped reference length: $unmapqrlength $unmapqrlengthper\%\n";
    print "Sequence IDs: ";
    scalar(@unmappedqr) < 1 ? print "None\n": print join(",",@unmappedqr),"\n";
    
    print "\nMapping details:\n";
    foreach my $key (@refseqID){
	my $coverpercent=sprintf("%.2f",$refcover{$key}/$reflength{$key}*100);
	my $idenpercent=sprintf("%.2f",$refidentity{$key});
	print "$key\t$reflength{$key}\t$refcover{$key}\t$coverpercent\%\t$idenpercent\%\t$refmapID{$key}\n";
    }
}
else{
    print "$queryfile\t",scalar(keys %querylength),"\t$querytotallength\t$totalcoverper\%\t$totalqueryper\%\t$totalsimper\%\t";
    print scalar(keys %refcover),"/",scalar(keys %reflength),"\t",$mappedqrnum,"/",scalar(keys %querylength),"\t","$unmapqrlength\t$unmapqrlengthper\%\n";
    my @outnums;
    foreach my $key (@refseqID){
	my $coverpercent=sprintf("%.2f",$refcover{$key}/$reflength{$key}*100);
	my $idenpercent=sprintf("%.2f",$refidentity{$key});
	push (@outnums,"$coverpercent\%\t$idenpercent");
    }
    print "$queryfile\t",join("\t",@outnums),"\n";
}
exit(1);