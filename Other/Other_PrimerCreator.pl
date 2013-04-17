#!/usr/bin/perl
use strict;

my $parFasfile=shift;
my $lastseqID;
my $seq;
my $primerlength=25;
open(my $fh_fasfile,"$parFasfile") || die "Can't open file $parFasfile\n";
print "Sequence_ID\tfwd_primer_seq\trev_primer_seq\tfwd_primer_name\trev_primer_name\n";
while(<$fh_fasfile>){
  chomp();
  if(/>(\S+)/){
    if($seq ne ""){
      my $pForwad=substr($seq,0,$primerlength);
      my $seqlength=length($seq);
      my $pEStart=$seqlength-$primerlength;
      my $pReverse=substr($seq,$pEStart);
      $pReverse=~tr/atgcATGC/tacgTACG/;
      $pReverse=reverse $pReverse;
      print "$lastseqID\t$pForwad\t$pReverse\t$lastseqID\_forward\t$lastseqID\_reverse\n";
    }
    $lastseqID=$1;
    $seq="";
  }
  else{
    $seq.=$_;
  }
}
close $fh_fasfile;
my $pForwad=substr($seq,0,$primerlength);
my $seqlength=length($seq);
my $pEStart=$seqlength-$primerlength;
my $pReverse=substr($seq,$pEStart);
$pReverse=~tr/atgcATGC/tacgTACG/;
$pReverse=reverse $pReverse;
print "$lastseqID\t$pForwad\t$pReverse\t$lastseqID\_forward\t$lastseqID\_reverse\n";

exit(1);
