#!/usr/bin/perl;
use strict;
use Getopt::Long;

use vars qw($parGTFFile $parExpFile $parOutDir $parGenomeInfo);
GetOptions("g=s" => \$parGTFFile, "e=s" => \$parExpFile, "m=s" => \$parGenomeInfo,"o=s" => \$parOutDir);
if($parGTFFile eq "" || $parExpFile eq "" || $parOutDir eq "" || $parGenomeInfo eq ""){
  print "Usage:$0 -g GTF_FILE -m GenomeInfo_FILE -e Expression_FILE -o Output_Directory\n";
  exit(0);
}

#my $windowSize=5000;
my $windowSize=100000;
my $frameSize=100000;
  
my %geneinfo;
my %geneExp;
my %chromGenes;
my %genomeExp;
my %genomeInfo;
open(my $fh_genomeinfo,$parGenomeInfo) || die "Can't open file: $parGenomeInfo\n";
while(<$fh_genomeinfo>){
  chomp();
  my @lines=split(/\t/,$_);
  $genomeInfo{$lines[0]}=$lines[1];
}
close $fh_genomeinfo;

&readGTFFile($parGTFFile);

open(my $fh_expfile,$parExpFile) || die "Can't open file: $parExpFile\n";
while(<$fh_expfile>){
  ##file format:    GeneID\tRPKM\n
  chomp();
  my @lines=split(/\t/,$_);
  $geneExp{$lines[0]}=$lines[1];
}
close $fh_expfile;

mkdir $parOutDir if(!(-e($parOutDir)));
open(my $outfileall,">$parOutDir/$parOutDir");

foreach my $chrom(keys %genomeInfo){
  next if(!exists($chromGenes{$chrom}));
  my $chrlength=$genomeInfo{$chrom};
  my $chromS=0;
  my $chromE=$chromS+$frameSize;
  my @chrGenes=@{$chromGenes{$chrom}};
  my $gstartIndex=0;
  open(my $outfile,">$parOutDir/$chrom.expr");
  while($chromE<$chrlength){
    my $totalExps=0;
    my $totalGenelength=0;
	my $isfirst=0;
    for(my $i=$gstartIndex;$i<@chrGenes;$i++){
      my $geneID=$chrGenes[$i];
      if(($geneinfo{$geneID}->{gstart}<=$chromS && $geneinfo{$geneID}->{gend} >=$chromS) ||($geneinfo{$geneID}->{gstart}>=$chromS && $geneinfo{$geneID}->{gend} <= $chromE) || ($geneinfo{$geneID}->{gstart}<=$chromE && $geneinfo{$geneID}->{gend} >= $chromE))
        {
          if(exists($geneExp{$geneID})){
            $totalExps+=$geneExp{$geneID};
            $totalGenelength+=$geneinfo{$geneID}->{glength};
            if($isfirst==0){
              $gstartIndex=$i;
              $isfirst=1;
            }
          }
        }
      else{
        last if($geneinfo{$geneID}->{gstart}>=$chromE);
      }
    }
    my $aveExp;
    if($totalExps!=0){
      $aveExp=$totalExps/$totalGenelength*1000;
    }
    my $aveExpOut=sprintf("%.4f",$aveExp);
    print $outfile "$chromS:$chromE\t$aveExpOut\t$totalExps\t$totalGenelength\n";
    my $regE=$chromS+$windowSize-1;
    my $lchrom="rn".lc($chrom);
    #my $logave=log($aveExpOut+0.01)/log(2);
    print $outfileall "$lchrom $chromS $regE $aveExpOut\n";
    $chromS+=$windowSize;
    $chromE=$chromS+$frameSize;
  }
  close $outfile;
}
close $outfileall;

exit(0);

sub readGTFFile{
  #this function only to record the exons of the genes;
  my $GTFFile=shift;
  my $lastgid;
  my $laststrand;
  my $lastgenetype;
  my $lastchrom;
  open(my $fh_gtffile,$GTFFile) || die "Can't open file: $GTFFile\n";
  my @exons;
  my @geneIDs;
  while(<$fh_gtffile>){
    my @lines=split(/\t/,$_);
    next if($lines[2] ne "exon");
    my $idinfos=$lines[8];
    my $geneid;
    if($lastchrom ne $lines[0]){
      if(scalar(@geneIDs)>0){
        my @chrgids=@geneIDs;
        $chromGenes{$lastchrom}=\@chrgids;
      }
      @geneIDs=();
      $lastchrom=$lines[0];
    }
    if($idinfos=~/gene_id \"(\S+)\";/){
      $geneid=$1;
    }
    if($geneid ne $lastgid){
      $geneinfo{$geneid}->{strand}=$lines[6];
      $geneinfo{$geneid}->{genetype}=$lines[1];
      $geneinfo{$geneid}->{chrom}=$lines[0];
      if(scalar(@exons)>0){
        my @exoninfo=&exonsort(@exons);
        $geneinfo{$lastgid}->{exons}=\@exoninfo;
        $geneinfo{$lastgid}->{exonnum}=scalar(@exoninfo);
        my($gs,$ge,$gl)=&genelength(@exoninfo);
        $geneinfo{$lastgid}->{gstart}=$gs;
        $geneinfo{$lastgid}->{gend}=$ge;
        $geneinfo{$lastgid}->{glength}=$gl;
        for(my $i=0;$i<@exoninfo;$i++){
          my $s=$exoninfo[$i]->{start};
          my $e=$exoninfo[$i]->{end};
        }
      }
      $lastgid=$geneid;
      push(@geneIDs, $geneid);
      @exons=();
      my %exoninfo;
      $exoninfo{start}=$lines[3];
      $exoninfo{end}=$lines[4];
      push(@exons,\%exoninfo);
    }
    else{
      my %exoninfo;
      $exoninfo{start}=$lines[3];
      $exoninfo{end}=$lines[4];
      push(@exons,\%exoninfo);
    }
  }
  {
    my @exoninfo=&exonsort(@exons);
    $geneinfo{$lastgid}->{exons}=\@exoninfo;
    $geneinfo{$lastgid}->{exonnum}=scalar(@exoninfo);
    my($gs,$ge,$gl)=&genelength(@exoninfo);
    $geneinfo{$lastgid}->{gstart}=$gs;
    $geneinfo{$lastgid}->{gend}=$ge;
    $geneinfo{$lastgid}->{glength}=$gl;
     my @chrgids=@geneIDs;
    $chromGenes{$lastchrom}=\@chrgids;
  }
  close $fh_gtffile;
}

sub genelength{
  my @exons=@_;
  my $glength;
  my $gs=$exons[0]->{start};
  my $ge;
  for(my $i=0;$i<@exons;$i++){
    my $es=$exons[$i]->{start};
    my $ee=$exons[$i]->{end};
    $ge=$ee;
    $glength+=$ee-$es+1;
  }
  return ($gs,$ge,$glength);
}

sub exonsort{
  my @exons=@_;
  my @sorted;
  my @fixed;
  while(scalar(@exons)>0){
    my $minindex;
    my $minstart;
    $minstart=$exons[0]->{start};
    for(my $i=1;$i<@exons;$i++){
      my $exonstart=$exons[$i]->{start};
      if($exonstart<$minstart){
        $minstart=$exonstart;
        $minindex=$i;
      }
    }
    push(@sorted,$exons[$minindex]);
    splice(@exons,$minindex,1);
  }
  
  for(my $i=0;$i<@sorted-1;$i++){
    my $reg1s=$sorted[$i]->{start};
    my $reg1e=$sorted[$i]->{end};
    my $reg2s=$sorted[$i+1]->{start};
    my $reg2e=$sorted[$i+1]->{end};
    
    if($reg1e>=$reg2s){
      if($reg1e<$reg2e){
        $sorted[$i]->{end}=$reg2e;
        splice(@sorted,$i+1,1);
        $i=-1;
      }
      else{
        splice(@sorted,$i+1,1);
      }
    }
  }  
  return @sorted;
}
