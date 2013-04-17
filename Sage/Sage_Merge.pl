#!/usr/bin/perl
use strict;

my ($confile,$outfile,$maxline)=@ARGV;
if(@ARGV<3){
    print "Usage:$0 Configure_File Output_File Maxline_ID\n";
    exit;
}

open(confile,$confile) || die "Can't open configure_file:$confile\n";
my @filenames;
my @colnames;
while(<confile>){
    chomp();
    push(@filenames,$_);
    my @fileN=split(/\./,$_);
    push(@colnames,$fileN[0]);
}
close;
my $matrix;

for(my $i=0;$i<@filenames;$i++){
    open(inputfile,$filenames[$i]);
    while(<inputfile>){
        chomp();
        my @lines=split(/\t/,$_);
        $matrix->[$lines[0]][$i]=$lines[1];
    }
    close inputfile;
}

open (outfile,">$outfile");
print outfile "$colnames[0]";
for (my $i=1;$i<@colnames;$i++){
    print outfile "\t$colnames[$i]";
}
print outfile "\n";

for(my $i=1;$i<=$maxline;$i++){
# print outfile "$i";
    for(my $j=0;$j<@colnames;$j++){
	if($j==0){
		$matrix->[$i][$j]>0 ? print outfile "$matrix->[$i][$j]" : print outfile "0";
	}
	else{
        	$matrix->[$i][$j]>0 ? print outfile "\t$matrix->[$i][$j]" : print outfile "\t0";
	}
    }
    print outfile "\n";
}
