#!/usr/bin/perl
#This file is used to add addtional information into the table according to the ID numbers;

use strict;
my ($inputfile,$idfield,$inforfile,$firstline)=@ARGV;
die "usage:$0 <input_file> <id_field> <add_info> [skip_first_line: default 1]\n" if(@ARGV<3);

$firstline=1 if($firstline eq "");
my %addinfo;
my $infonumber;
my $title;

open(my $fh_inforfile,$inforfile) || die "Can't open file:$inforfile\n";
my $headline=<$fh_inforfile>;
chomp($headline);
my @titles=split(/\t/,$headline);
shift @titles;
$title=join("\t",@titles);
while(<$fh_inforfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my $id=shift(@lines);
    $infonumber=scalar(@lines);
    $addinfo{$id}=join("\t",@lines);
}
close $fh_inforfile;
my $undef="undef";
for(my $i=1;$i<$infonumber;$i++){
    $undef.="\tundef";
}

open(my $fh_inputfile,$inputfile) || die "Can't open file:$inputfile\n";
if($firstline==1){
    my @heads=split(/\t/,<$fh_inputfile>);
    splice(@heads,$idfield,0,$title);
    print join("\t",@heads);
}

while(<$fh_inputfile>){
    chomp();
    $_=~s/\"//g;
    my @lines=split(/\t/,$_);
    if(exists($addinfo{$lines[$idfield]})){
        splice(@lines,$idfield+1,0,$addinfo{$lines[$idfield]});
        print join("\t",@lines),"\n";
    }
    else{
        splice(@lines,$idfield+1,0,$undef);
        print join("\t",@lines),"\n";
    }
}
close $fh_inputfile;

exit(1);
