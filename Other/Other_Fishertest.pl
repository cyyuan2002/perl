#!/usr/bin/perl
use strict;
my ($input11,$input12,$input21,$input22)=@ARGV;
if(@ARGV<4){
	print "Usage:$0 A11 A12 A21 A22\n";
	exit(0);
}
my $pvalue=&fishertest($input11,$input12,$input21,$input22);
print "$input11 $input12\n$input21 $input22\np-value:$pvalue\n";


sub fishertest{
  my ($n11,$n12,$n21,$n22)=@_;
  use Text::NSP::Measures::2D::Fisher2::right;
  my $n1p=$n11+$n12;
  my $np1=$n11+$n21;
  my $npp=$n11+$n12+$n21+$n22;
  my $right_value=calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
  return $right_value;
}
