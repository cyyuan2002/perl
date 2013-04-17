#!/usr/bin/perl
use strict;

my ($parGTFFile,$parOutFile)=@ARGV;
if(@ARGV<2){
  print "Usage:$0 GTF_File Output_File\n";
  exit(0);
}

my %geneinfo;
my %chromGenes;

&readGTFFile($parGTFFile);
open(my $fh_outfile,">$parOutFile");

foreach my $chrom(keys %chromGenes){
  my @chrGenes=@{$chromGenes{$chrom}};
  for(my $i=0;$i<@chrGenes;$i++){
    my %geneinfo=%{$geneinfo{$chrGenes[$i]}};
    print $fh_outfile "$chrom\t$geneinfo{gstart}\t$geneinfo{gend}\t$geneinfo{strand}\t$chrGenes[$i]\t$geneinfo{glength}\t$geneinfo{exonnum}\t";
    my @exons=@{$geneinfo{exons}};
    my $exonsinfo;
    foreach my $exoninfo(@exons){
      my $es=$exoninfo->{start};
      my $ee=$exoninfo->{end};
      $exonsinfo.="$es:$ee,";
    }
    print $fh_outfile "$exonsinfo\n";
  }
}
close $fh_outfile;

sub readGTFFile{
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
