#!/usr/bin/perl
use strict;

my ($parFile1,$parField1,$parFile2,$parField2,$parOutFile) = @ARGV;
if(@ARGV<5){
  print "Usage:$0 IDFile IDField DataFile IDField OutputFile\n";
  exit(0);
}

my %ids;
open(my $fh_file1,$parFile1) || die "Can't open file: $parFile1\n";
while(<$fh_file1>){
  chomp();
  my @lines=split(/\t/,$_);
  $ids{$lines[$parField1]}=1;
}
close $fh_file1;

open(my $fh_file2,$parFile2) || die "Can't open file: $parFile2\n";
open (my $outfile,">$parOutFile");
while(<$fh_file2>){
  chomp();
  my @lines=split(/\t/,$_);
  my $idN=$lines[$parField2];
  $idN=~s/\"//g;
  if(exists($ids{$idN})){
    print $outfile "$_\n";
  }
}
close $fh_file2;

exit(1);
