#!/usr/bin/perl

use 5.006;
use strict;
use warnings;
use Getopt::Long;

=head1 Name

GenRepeat.pl -- simulate repeats on linux system.

=head1 Description

This program is used to simulate repeats.
There are three kinds of repeats in genome:
(1)simple tandem repeats,STRs <300bp 
(2)variable number of tandem repeat,VNTR >300bp 
(3)Transposon	>500bp

=head1 Version

  Author: Zhan dongliang, zhandongliang@genomics.org.cn
  Version: 1.0,  Date: 2010-11-22

=head1 Usage
  
  perl GenRepeat.pl
  --size <num>      the genome size of 	DNA[100000000]
  --STR_node <num>  set number of types in simple tadem repeats[0]
  --STR_len  <num>  set the length of simple tadem repeats in the sequece[STR_node * 300]
  --VNTR_node <num> set number of types in minisatelite[0]
  --VNTR_len  <num> set the length of minisatelite in the sequece[VNTR_node *1000]
  --TRAN_node <num> set number of types in transposon[0]
  --TRAN_min <num>  set minimal size of transposons [500]
  --TRAN_max <num>  set minimal size of transposons [1000]
  --TRAN_len  <num> set the length of transposons in the sequece[TRAN_len * TRAN_node]
  --Mode <num>      set the mode of generation of repeats[1]
                    1:insert sequence between two repeats
                    2:replace part of sequence by repeats
                    3:replace part of sequence by repeats,but won't replace the repeat,only for repeats < 5%
  --help            output help information to screen  

=head1 Exmple

  perl GenRepeat.pl --size 100000000 --STR_node 1000 --STR_len 300000 --VNTR_node 500 --VNTR_len 500000 
                    --TRAN_node 100 --TRAN_min 500 --TRAN_max 1000 --TRAN_len 100000 --Mode 1

=cut
my($size,$str_node,$str_len,$vntr_node,$vntr_len,$tran_node,$tran_min,$tran_max,$tran_len,$mode,$help);
die `pod2text $0` if (@ARGV == 0);
GetOptions(
	"size:i"=>\$size,
	"STR_node:i"=>\$str_node,
	"STR_len:i"=>\$str_len,
	"VNTR_node:i"=>\$vntr_node,
	"VNTR_len:i"=>\$vntr_len,
	"TRAN_node:i"=>\$tran_node,
	"TRAN_min:i"=>\$tran_min,
	"TRAN_max:i"=>\$tran_max,
	"TRAN_len:i"=>\$tran_len,
	"Mode:i"=>\$mode,
	"help"=>\$help
);
die `pod2text $0` if ($help);

#initial parameter
$mode ||= 1;
$str_node ||= 0;
$str_len ||= 300*$str_node;
$vntr_node ||=0;
$vntr_len ||= 1000*$vntr_node;
$tran_node ||=0;
$tran_min ||= 500;
$tran_max ||= 1000;
$tran_len ||= $tran_max*$tran_node;
$size ||= 100000000;

my $gap=0;
my $t_len=0;
#need to be modify if your want to generate more repeats
sub replace{
	our @LIST;
	my ($ref,$start,$seq) = @_;
	my $len = length($seq);
	my $size = length(${$ref});
	my $end = $start+$len-1;
	if($start+$len>=$size){
		return -1;
	}
	foreach(@LIST){
		my @array = split(/\t/,$_);
		if($start>=$array[0]&&$start<=$array[1]||$end>=$array[0]&&$end<=$array[1]){
			return -1;
		}
	}
	my $info = "$start\t$end";
	push(@LIST,$info);
	substr(${$ref},$start,$len) = $seq;
	return 1;
}

#generate random sequence 
sub genRandomSeq{
	my $seq_num = shift;
	my $seq;
	my $count = 0;
      	while($count!=$seq_num){
        	my $ran=rand;
         	$ran = $ran*10000%4;
        	my $key = '';
         	if($ran ==0){
              		$key = 'A';
         	}elsif($ran ==1){
              		$key ='T';
         	}elsif($ran ==2){
              		$key ='G'
         	}else{
              		$key ='C';
        	}
         	$count++;
		$seq .= $key;
      	}
	return $seq;
}

#generate simle tandem repeats  len:300*num
sub genSimpleTandem{
	my($num,$len) = @_;
	if($num==0||$len==0){
		return;
	}
	open REPEAT,">>lardo.repeat";
	my $avg = $len/$num;
	my $total = 0;
	for(my $i=0;$i<$num;$i++){
		#max len=7*40 280
		my $seed = int(2+rand(6));
		my $repeat = genRandomSeq($seed);
		my $t = int(2+rand(40));
		$repeat = $repeat x$t;
		my $leng = length($repeat);
		my $n = int(($avg/$leng)+0.5);
		$gap += $n;
		print REPEAT "$repeat\t$n\n";
		$total += $leng*$n;
	}
	close REPEAT;
	$t_len += $total;
	print STDERR "SIMPLE TANDOM REPEAT:$total\n";
}

#generate mini tandem repeats len:1000*num
sub genMiliTandem{
	my($num,$len) = @_;
	if($num==0||$len==0){
		return;
	}
	open REPEAT,">>lardo.repeat";
	my $avg = $len/$num;
	my $total = 0;
	for(my $i=0;$i<$num;$i++){
		#max len=1000
		my $seed = int(10+rand(15));
		my $repeat = genRandomSeq($seed);
		my $t = int(2+rand(40));
		$repeat = $repeat x$t;
		my $leng = length($repeat);
		my $n = int(($avg/$leng)+0.5);
		$gap += $n;
		print REPEAT "$repeat\t$n\n";
		$total += $leng*$n;
	}
	close REPEAT;
	$t_len += $total;
	print STDERR "MILI TANDOM REPEAT:$total\n";
}

#generate transporter len:$max*num
sub genTransport{
	my($num,$len,$min,$max) = @_;
	if($num==0||$len==0){
		return;
	}
	open REPEAT,">>lardo.repeat";
	my $avg = $len/$num;
	my $total = 0;
	for(my $i=0;$i<$num;$i++){
		my $seed = int($min+rand($max-$min));
		my $repeat = genRandomSeq($seed);
		my $leng = length($repeat);
		my $n = int(($avg/$leng)+0.5);
		$gap += $n;
		print REPEAT "$repeat\t$n\n";
		$total += $leng*$n;
	}
	close REPEAT;
	$t_len += $total;
	print STDERR "Transport REPEAT:$total\n";
}


##################### main ####################################################

genTransport($tran_node,$tran_len,$tran_min,$tran_max);
genMiliTandem($vntr_node,$vntr_len);
genSimpleTandem($str_node,$str_len);
my $percent = $t_len/$size*100;
#information about repeats
print STDERR "num:$gap\n";
print STDERR "repeat:$t_len\n";
print STDERR "REPEAT:$percent%\n";
print STDERR "mode:$mode\n";

my $avg = int (($size-$t_len)/($gap+1));
print STDERR "average sequence len:$avg\n";


open OUT,">lardo.fa";
print OUT ">LARDO\n";
unless(-e "lardo.repeat"){
	my $ref = genRandomSeq($size);
	print OUT "$ref\n";
}
open IN,"lardo.repeat" or die "finish\n";
my $count = 1;
my $ref = "";
if($mode!=1){
	$ref = genRandomSeq($size);
}
while(<IN>){
	chomp;
	my @array = split(/\t/,$_);
	my $seq = $array[0];
	my $num = $array[1];
	for(my $i=0;$i<$num;$i++){
		if($mode eq 3){
			my $start = int rand($size);
			if(replace(\$ref,$start,$seq)==-1){
				redo;
			}
		}elsif($mode eq 1){
			my $tem = genRandomSeq($avg);
			$ref = "";
			$ref .= $tem;
			$ref .= $seq;
			print OUT "$ref";
		}elsif($mode eq 2){
			my $start = int rand($size);
			my $len = length($seq);
			substr($ref,$start,$len) = $seq;
		}
	}
=f
	if($count%100 == 0){
		print STDERR "$count\n";
	}
=cut
	$count++;
}
if($mode!=1){
	print OUT "$ref";
}
close OUT;
close IN;
print STDERR "finished\n";

