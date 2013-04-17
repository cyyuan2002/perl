#! /usr/bin/perl
use strict;
#use warnings;


die "Version: 1.1\nData: 2012-05-08\nUsage:$0 <lane.lst> <libsize> [Q_Shift:default 64] [maxjob, default 10]\n" if @ARGV<2;

my $filter_gz_folder="~/bin/Filter_data_gz";
my $lane_lst = shift;
my $lib_size = shift;
my $Q_shift = shift;
my $Max_job = shift;

$Q_shift ||=64;
$Max_job ||=10;


my $file_name = `basename $lane_lst`;
chomp $file_name;

open OUT1, ">$file_name.filter.sh" or die "$!";
open OUT2, ">$file_name.dup.sh" or die "$!";
open OUT3, ">$file_name.stat.sh" or die "$!";


if ( -e "$file_name.stat.xls" ){
	my $tag = time;
	print STDERR "$file_name.stat.xls  already exists! \n";
	print STDERR "mv  $file_name.stat.xls $file_name.stat.xls.$tag \n";
	`mv  $file_name.stat.xls $file_name.stat.xls.$tag`;
}

open IN,$lane_lst or die "$!";

my $stat_output = '';
my $num = 0;
while(<IN>){
	chomp;
	$num++;
	
	my ($f1,$start1,$end1,$B_cutoff) = split(/\s+/);
	$start1 ||= 0;
	$end1 ||= 0;
	$B_cutoff ||= 40;
	
	my $line2 = <IN>;
	chomp $line2;
	my ($f2,$start2,$end2,$N_num) = split(/\s+/,$line2);
	$start2 ||= 0;
	$end2 ||= 0;
	if(not defined $N_num or $N_num eq "")
	{$N_num ||= 10;}
	
	my $name = `basename $f1`;
	chomp $name;
	my $name2 = `basename $f2`;
	chomp $name2;

	print OUT1 "$filter_gz_folder/filter_data_gz -y -z -Q $Q_shift -w $N_num -B $B_cutoff -l $lib_size -a $start1 -b $end1 -c $start2 -d $end2 $f1 $f2 $name.reads.stat $name.clean.tmp $name2.clean.tmp && mv $name.clean.tmp $name.clean && mv $name2.clean.tmp $name2.clean";

	print OUT2 "$filter_gz_folder/duplication $name.clean $name2.clean $name.clean.dup.clean $name2.clean.dup.clean $name.clean.dup.stat";
	
	my $parameter="$start1"."_"."$end1"."_"."$start2"."_"."$end2"."_"."$B_cutoff"."_"."$N_num";
	
	print OUT3 "$filter_gz_folder/stat.pl $lib_size $name.reads.stat $name.clean.dup.stat $parameter >>$file_name.stat.xls";
	

}
close IN;

print STDERR "Runing low quality filtering ... \n";
print STDERR "Finished. \n";
print STDERR "Runing duplicate filtering ... \n";
print STDERR "Finished. \n";
print STDERR "Runing stat. ... \n";
print STDERR "All finished. \n";
