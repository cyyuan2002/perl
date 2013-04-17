#!/usr/bin/perl
use strict;
use Compress::Zlib;

my $file=shift;
my $fh=gzopen($file,"rb");
my $i;
my $NCount;
my $linecount;
while($fh->gzreadline(my $line)){
    $i++;
    if($i==2){
	my $base=substr($line,0,1);
	if($base eq "N"){
	    $NCount++;
	    #print "$base\n";
	}
	$linecount++;
    }
    elsif($i==4){
	$i=0;
    }
}
$fh->gzclose();
my $percent=($NCount/$linecount)*100;
print "Total N: $NCount ; percentage: $percent%\n";