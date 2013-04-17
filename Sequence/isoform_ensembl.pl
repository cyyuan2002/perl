#!/usr/bin/perl

use strict;

my $file=shift;

open(file,$file) || die "Can't open file $file\n";
my $seqname;
my $seq;
my $outfile=$file.".filteredid";


open (outfile,">$outfile");

my %tids;
my %types;
my %length;
my %lines;

while(<file>){
	chomp();
	if(/^>/){
		if($seqname ne ""){
                    
		    my $length=length($seq);
                    $seqname=~s/>//g;
                    my @seqN=split(/\s/,$seqname);
                    my @type=split(/\:/,$seqN[1]);
                    my @gname=split(/\:/,$seqN[3]);
                    if($tids{$gname[1]} eq ""){
                        $tids{$gname[1]}=$seqN[0];
                        $types{$gname[1]}=$type[1];
                        $length{$gname[1]}=$length;
                        $lines{$gname[1]}=$seqname;
                    }
                    else{
                        if($types{$gname[1]} ne "known" && $types{$gname[1]} eq "known"){
                            $tids{$gname[1]}=$seqN[0];
                            $types{$gname[1]}=$type[1];
                            $length{$gname[1]}=$length;
                            $lines{$gname[1]}=$seqname;
                        }
                        else{
                            if($length{$gname[1]}<$length){
                                $tids{$gname[1]}=$seqN[0];
                                $types{$gname[1]}=$type[1];
                                $length{$gname[1]}=$length;
                                $lines{$gname[1]}=$seqname;
                            }
                        }
                    }
		}
		$seqname=$_;
		$seq="";
	}
	else{
		$seq=$seq.$_;
	}
}
close file;
my $length=length($seq);
$seqname=~s/>//g;
my @seqN=split(/\s/,$seqname);
my @type=split(/\:/,$seqN[1]);
my @gname=split(/\:/,$seqN[3]);
if($tids{$gname[1]} eq ""){
    $tids{$gname[1]}=$seqN[0];
    $types{$gname[1]}=$type[1];
    $length{$gname[1]}=$length;
    $lines{$gname[1]}=$seqname;
}
else{
    if($types{$gname[1]} ne "known" && $types{$gname[1]} eq "known"){
        $tids{$gname[1]}=$seqN[0];
        $types{$gname[1]}=$type[1];
        $length{$gname[1]}=$length;
        $lines{$gname[1]}=$seqname;
    }
    else{
        if($length{$gname[1]}<$length){
            $tids{$gname[1]}=$seqN[0];
            $types{$gname[1]}=$type[1];
            $length{$gname[1]}=$length;
            $lines{$gname[1]}=$seqname;
        }
    }
}

foreach my $key(keys %tids){
    print outfile "$lines{$key}\n";
}
close outfile;
