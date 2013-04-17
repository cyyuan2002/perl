##This program is used to remove the redundant lines in the 2 column table, merge the same id in column 1, and remove lines with same id in column2
#!/usr/bin/perl
use strict;

my $parFilein=shift;
my %id1;
my %id2;
open(my $fh_filein,$parFilein) || die "Can't open file $parFilein\n";
while(<$fh_filein>){
  chomp();
  my @lines=split(/\t/,$_);
  if(!exists($id1{$lines[0]})){
    if(!exists($id2{$lines[1]})){
      $id1{$lines[0]}=$lines[1];
      $id2{$lines[1]}=$lines[0];
    }
    else{
      delete($id1{$id2{$lines[1]}}) if(exists($id1{$id2{$lines[1]}}));
    }
  }
  else{
    next if($id1{$lines[0]} eq $lines[1]);
    $id1{$lines[0]}.=",$lines[1]";
  }
}
close $fh_filein;

foreach my $key(keys %id1){
  print "$key\t$id1{$key}\n";
}
