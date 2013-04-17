#!/usr/bin/perl
use strict;

my ($file1,$file2,$outfile)=@ARGV;
my @fileinfor1;
my @fileinfor2;
my @orthologs;
my @searched;

open (file1,$file1) || die "Can't open file: $file1\n";
@fileinfor1=<file1>;
close file1;

open (file2,$file2) || die "Can't open file: $file2\n";
@fileinfor2=<file2>;
close file2;

for(my $i=0;$i<@fileinfor1;$i++){
    my @lines1=split(/\t/,$fileinfor1[$i]);
    my $queryN=$lines1[0];
    my $subN=$lines1[1];
    my $isskip=0;
    for(my $k=0;$k<@searched;$k++){
        if($queryN eq $searched[$k]){
            $isskip=1;
            last;
        }
    }
        
    if ($isskip==1){
        next;
    }
    
    for (my $j=0;$j<@fileinfor2;$j++){
        my @lines2=split(/\t/,$fileinfor2[$j]);
        if($lines2[0] eq $subN){
            if($lines2[1] eq $queryN){
                my $ortholog=$queryN."\t".$subN;
                push(@orthologs,$ortholog);
                push (@searched,$subN);
            }
            else{
                last;
            }
        }
    }
}

open(out,">$outfile");
for (my $i=0;$i<@orthologs;$i++){
    print out "$orthologs[$i]\n";
}
