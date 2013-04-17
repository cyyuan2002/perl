#!/usr/bin/perl
use strict;

my ($parExpFile,$parAnnoFile,$parGeneList,$parToppercent,$parOutFile)=@ARGV;
if(@ARGV<5){
  print "usage:$0 SortedExpressFile AnnotationFile GeneList Toppercent OutputFile\n";
  exit(0);
}

my %geneinfo;
my %chromgene;
my $windowsize=100000;
open(my $fh_annofile,"$parAnnoFile") || die "Can't open file: $parAnnoFile\n";
my $lastchrom;
my @geneIDs;
while(<$fh_annofile>){
  chomp();
  my @lines=split(/\t/,$_);
  $geneinfo{$lines[4]}->{gstart}=$lines[1];
  $geneinfo{$lines[4]}->{gend}=$lines[2];
  $geneinfo{$lines[4]}->{chrom}=$lines[0];
  if($lastchrom ne $lines[0]){
    if(scalar(@geneIDs)>0){
      my @geneNs=@geneIDs;
      $chromgene{$lastchrom}=\@geneNs;
    }
    $lastchrom=$lines[0];
    @geneIDs=();
  }
  push(@geneIDs,$lines[4]);
}
my @geneNs=@geneIDs;
$chromgene{$lastchrom}=\@geneNs;
close $fh_annofile;

my %topgenes;
my @chroms=("2L","2R","3L","3R","4");

open(my $fh_expfile,"$parExpFile") || die "Can't open file: $parExpFile\n";
my @expInfos=<$fh_expfile>;
close $fh_expfile;
my $topNum=int(scalar(@expInfos)*$parToppercent/100);
for(my $i=0;$i<$topNum;$i++){
  my $lineinfo=$expInfos[$i];
  my @lines=split(/\s/,$lineinfo);
  my $chrom=&transchrom($lines[0]);
  #print "$lines[0]\t$chrom\n";
  next if($chrom eq "na");
  my $chromS=$lines[1];
  my $chromE=$chromS+$windowsize;
  my @chrGenes=@{$chromgene{$chrom}};
  my $geneNum=scalar(@chrGenes);
  for(my $i=0;$i<@chrGenes;$i++){
      my $geneID=$chrGenes[$i];
     # my $gs=$geneinfo{$geneID}->{gstart};
     # my $ge=$geneinfo{$geneID}->{gend};
      if(($geneinfo{$geneID}->{gstart}<=$chromS && $geneinfo{$geneID}->{gend} >=$chromS) ||($geneinfo{$geneID}->{gstart}>=$chromS && $geneinfo{$geneID}->{gend} <= $chromE) || ($geneinfo{$geneID}->{gstart}<=$chromE && $geneinfo{$geneID}->{gend} >= $chromE))
        {
          $topgenes{$geneID}=1;
        }
      else{
        last if($geneinfo{$geneID}->{gstart}>=$chromE);
      }
    }
}

open(my $fh_genelist,"$parGeneList");
my $overcount=0;
my $genelistcount=0;
while(<$fh_genelist>){
  chomp();
  $genelistcount++;
  my @lines=split(/\t/,$_);
  if(exists($topgenes{$lines[0]})){
    $overcount++;
  }
}
close $fh_genelist;

open(my $fh_fileout,">$parOutFile");
foreach my $gid(keys %topgenes){
  my %geneI=%{$geneinfo{$gid}};
  print $fh_fileout "$gid\t$geneI{chrom}\n";
}
close $fh_fileout;

my @topgenes=keys(%topgenes);
my $topgenenum=scalar(@topgenes);

my $totalgenenum;
foreach my $chrom(@chroms){
  $totalgenenum+=scalar(@{$chromgene{$chrom}});
}
print "Top $parToppercent\% genes: $topgenenum\n";
print "Overlapped genes: $overcount\n";
my $pvalue=&fishertest($overcount,$genelistcount,$topgenenum,$totalgenenum);
print "Total genes: $totalgenenum\n";
print "Pvalue: $pvalue\n";

sub transchrom{
  my $chrom=shift;
  if($chrom eq "rn2l"){
    return "2L";
  }
  elsif($chrom eq "rn2r"){
    return "2R";
  }
  elsif($chrom eq "rn3l"){
    return "3L";
  }
  elsif($chrom eq "rn3r"){
    return "3R";
  }
  elsif($chrom eq "rn4"){
    return "4";
  }
  elsif($chrom eq "rnx"){
    return "X";
  }
  else{
    return "na";
  }
}

sub fishertest{
  my ($n11,$n12,$n21,$n22)=@_;
  use Text::NSP::Measures::2D::Fisher2::right;
  my $n1p=$n11+$n12;
  my $np1=$n11+$n21;
  my $npp=$n11+$n12+$n21+$n22;
  my $right_value=calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
  return $right_value;
}
