#!/usr/bin/perl

#===============================================================================
#
#         FILE: 
#
#        USAGE:
#
#  DESCRIPTION:This scipt is used to offset the mapping quality of Tophat from 255 to -1
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 
#     REVISION:
#===============================================================================

use strict;

my $BamFile=shift;
my $samtools=`which samtools`;
if($samtools eq ""){
    print "Can't find samtools\n";
    exit(-1);
}

`samtools view -h -o $BamFile.sam $BamFile`;
open(my $fh_samfile, "<", "$BamFile.sam") or die "Can't find file $BamFile.sam\n";
open(my $fh_outsam, ">", "$BamFile.edt.sam");
while (<$fh_samfile>) {
    chomp();
    my @lines=split(/\t/);
    if ($lines[5] ==255) {
        $lines[5]=-1;
        print $fh_outsam join("\t",@lines);
    }
    else{
        print $fh_outsam $_;
    }
    print $fh_outsam "\n";
}
close $fh_samfile;

`samtools view -Sb $BamFile.edt.sam > $BamFile.offset.bam`;

unlink "$BamFile.sam";
unlink "$BamFile.edt.sam";
