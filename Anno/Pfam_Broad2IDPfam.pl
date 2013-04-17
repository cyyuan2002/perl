#!/usr/bin/perl
##This program is used to deal with pfam file from broad
use strict;

my ($file,$idfile)=@ARGV;
die "usage:$0 <pfam_file> <geneids>\n" if(@ARGV<2);


open(my $fh_file,$file) || die "Can't open file: $file\n";

my $lastlocus;
my $proname;
my %pfamname;
my %pfamdescript;

my %pronames;
my %pfamnames;
my %pfamdescripts;

<$fh_file>;

while(<$fh_file>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastlocus ne $lines[1]){
        if($lastlocus ne ""){
            $pronames{$lastlocus}=$proname;
            $pfamnames{$lastlocus}=join(" | ",(keys %pfamname));
            $pfamdescripts{$lastlocus}=join(" | ",(keys %pfamdescript));
        }
        $lastlocus=$lines[1];
        $proname=$lines[0];
        %pfamname=();
        $pfamname{$lines[4]}=1;
        %pfamdescript=();
        $pfamdescript{$lines[5]}=1;
    }
    else{
        if(!exists($pfamname{$lines[4]})){
            $pfamname{$lines[4]}=1;
            $pfamdescript{$lines[5]}=1;
        }
    }
}
close $fh_file;
$pronames{$lastlocus}=$proname;
$pfamnames{$lastlocus}=join("|",(keys %pfamname));
$pfamdescripts{$lastlocus}=join("|",(keys %pfamdescript));

print "LOCUS\tprotein_name\tPfam_name\tPfam_desc\n";
my @ids;
open(my $fh_idfile,$idfile) || die "Can't open file:$idfile\n";
while(<$fh_idfile>){
    chomp();
    push(@ids,$_);
}
close $fh_idfile;

foreach my $key(@ids){
    if(exists($pronames{$key})){
        print "$key\t$pronames{$key}\t$pfamnames{$key}\t$pfamdescripts{$key}\n";
    }
    else{
        print "$key\tundef\tundef\tundef\n";
    }
}

exit(1);
