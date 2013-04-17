#!/usr/bin/perl
use strict;

my ($parInputFiles,$parOutputFile)=@ARGV;
if(@ARGV<2){
  print "Usage:$0 GeneExp_File(example: fileA,fileB,...) OutputFile\n";
  exit(0)
}

my @files=split(/,/,$parInputFiles);

foreach my $fileN(@files){
  if(!-e($fileN)){
    print "Can't find file: $fileN\n";
    exit(0);
  }
}

my %expinfo;
foreach my $fileN(@files){
  open(my $fh_file,$fileN);
  <fh_file>;
  while(<$fh_file>){
    chomp();
    my @lines=split(/\t/,$_);
    my $geneid=$lines[0];
    my $expvalue=$lines[9];
    if(exists($expinfo{$geneid})){
      $expinfo{$geneid}+=$expvalue;
    }
    else{
      $expinfo{$geneid}=$expvalue;
   } 
  }
}

my $sampleNum=scalar(@files);
open(my $fh_outfile,">$parOutputFile");
foreach my $key(keys %expinfo){
  my $expave=$expinfo{$key}/$sampleNum;
  my $expvalout=sprintf("%.4f",$expave);
  print $fh_outfile "$key\t$expvalout\n";
}

close $fh_outfile;
