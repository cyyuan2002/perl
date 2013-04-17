#!/usr/bin/perl
use strict;

my ($parFile1,$parField1,$parFile2,$parField2,$parFieldMatch,$parIsRed)=@ARGV;
if(@ARGV<5){
  print "Usage:$0 InputFile Field2Convert TableFile Field2Match Field2Output IsRed(0:no,1:yes)\n";
  exit(0);
}
$parIsRed=0 if($parIsRed eq "");

my %idNames;
my %idConv;
open(my $fh_file1,$parFile1) || die "Can't open file: $parFile1\n";
while(<$fh_file1>){
  chomp();
  my @lines=split(/\t/,$_);
  $idNames{$lines[$parField1]}=1;
}
close $fh_file1;

open(my $fh_file2,$parFile2) || die "Can't open file: $parFile2\n";
while(<$fh_file2>){
  chomp();
  my @lines=split(/\t/,$_);
  if(exists($idNames{$lines[$parField2]})){
    if(exists($idConv{$lines[$parField2]})){
      if($parIsRed==1){
        $idConv{$lines[$parField2]}.=",$lines[$parFieldMatch]";
      }
    }
    else{
      $idConv{$lines[$parField2]}=$lines[$parFieldMatch];
    }
  }
}
close $fh_file2;

foreach my $key(keys %idConv){
  print "$idConv{$key}\t$key\n";
}

exit(1);
