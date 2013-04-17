#!/usr/bin/perl
use strict;

my ($VCFfile,$SV_TYPE,$CountFile,$LengthFile)=@ARGV;
open(my $fh_vcffile,"$VCFfile") || die "$!\n";

my $fh_countfile;
my $fh_lengthfile;

if(-e $CountFile){
    open($fh_countfile,">>$CountFile");
}
else{
    open($fh_countfile,">$CountFile");
}

if(-e $LengthFile){
    open($fh_lengthfile, ">>$LengthFile");
}
else{
    open($fh_lengthfile,">$LengthFile");
}
print $fh_lengthfile "$VCFfile";

my $SVCount=0;
my $TotalLength=0;
while(<$fh_vcffile>){
    next if(/^#/);
    chomp();
    $SVCount++;
    my @lines=split(/\t/,$_);
    my ($SV_s,$SV_e,$SV_length);
    if($SV_TYPE ne "CTX" && $SV_TYPE ne "ITX"){
        $SV_s=$lines[1];
        $lines[7]=~/END=(\d+)\;/;
        $SV_e=$1;
        $SV_length=abs($SV_e-$SV_s);
        print $fh_lengthfile "\t$SV_length";
        $TotalLength+=$SV_length;
    }
}
close $fh_vcffile;
if($SV_TYPE ne "CTX" && $SV_TYPE ne "ITX"){
    print $fh_lengthfile "\n";
    my $averlength=$TotalLength/$SVCount;
    print $fh_countfile "$VCFfile\t$SVCount\t$averlength\n";
}
else{
    print $fh_countfile "$VCFfile\t$SVCount\n";
}

