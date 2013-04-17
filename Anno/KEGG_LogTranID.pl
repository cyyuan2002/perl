#!/usr/bin/perl
use strict;

my ($parInfile,$parIdfile,$parFieldNum)=@ARGV;
if(@ARGV<2){
  print "Usage:$0 KEGG_LogFile Id_File Field_Num\n";
  exit(0);
}

my %ids;

open(my $fh_idfile,$parIdfile) || die "Can't open file: $parIdfile\n";
while(<$fh_idfile>){
  chomp();
  my @lines=split(/\t/,$_);
  $ids{$lines[0]}=$lines[1];
}
close $fh_idfile;

open(my $fh_infile,$parInfile) || die "Can't open file: $parInfile\n";
while(<$fh_infile>){
  chomp();
  my @lines=split(/\t/,$_);
  my @gnames=split(/,/,$lines[$parFieldNum]);
  for(my $i=0;$i<@gnames;$i++){
    $gnames[$i]=$ids{$gnames[$i]};
  }
  my $newids=join(",",@gnames);
  $lines[$parFieldNum]=$newids;
  print join("\t",@lines),"\n";
}
close $fh_infile;

exit(1);
