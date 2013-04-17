#!/usr/bin/perl
#this program is used map the tag information to the statistical results

use strict;

my ($mapfile,$comfile,$outfile) = @ARGV;
if(@ARGV<2){
    print stderr "usage:$0 Map_File Comp_File [Out_File]\n";
    exit(0);
}
if(!(-e $mapfile)){
    print stderr "Can't find file $mapfile\n";
    exit(0);
}
if(!(-e $comfile)){
    print stderr "Can't find file $comfile\n";
    exit(0);
}

if ($outfile eq ""){
    $outfile="$comfile.map";
}

my %Taginfo;
my %Taggene;
open (mapfile,$mapfile);
while(<mapfile>){
    chomp();
    my @lines=split(/\t/,$_);
    $Taginfo{$lines[0]}="$lines[2]\t$lines[3]\t$lines[4]";
    $Taggene{$lines[0]}=$lines[1];
}
close mapfile;

open (comfile,$comfile);
open (outfile,">$outfile");
while(<comfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my $tag=$lines[0];
    $tag=~s/\"//g;
    if(exists($Taginfo{$tag})){
	print outfile "$tag\t$Taggene{$tag}";
	for(my $i=1;$i<@lines;$i++){
	    print outfile "\t$lines[$i]";
	}
	print outfile "\t$Taginfo{$tag}\n";
    }
    else{
	print outfile "$tag\tNA";
	for(my $i=1;$i<@lines;$i++){
	    print outfile "\t$lines[$i]";
	}
	print outfile "\n";
    }
}
close comfile;
close outfile;
