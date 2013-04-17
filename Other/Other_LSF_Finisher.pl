#!/usr/bin/perl
#===============================================================================
#
#         FILE: Other_LSF_Finisher.pl
#
#        USAGE:
#
#  DESCRIPTION: This script is used to check the finished jobs on Broad_LSF
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION:
#      CREATED: 12/25/2012
#     REVISION:
#===============================================================================

use strict;
my $Suffix=shift;
$Suffix ||="log";

my @files=<*.$Suffix>;
my $outfile="Job_Finisher.sh";
my $isoutput=0;
open (my $fh_outfile,">$outfile");

my $count=0;
foreach my $file (@files){
    open(my $fh_file,"$file");
    $count++;
    my $iskilled=0;
    my $command="";
    while(<$fh_file>){
        if(/^Subject: Job \d+: \<(.+)\> Exited/){
            $isoutput=1;
            $iskilled=1;
            $command=$1;
            print stderr "\n$file\n";
            print $fh_outfile "rm $file\n";
        }
        elsif(/^TERM_RUNLIMIT:/){
            print $fh_outfile "bsub -o $file -q week \"$command\"\n";
            print stderr "RUNLIMIT: $command\n";
        }
        elsif(/^TERM_MEMLIMIT:/){
            print $fh_outfile "bsub -o $file -R \"rusage\[mem=4\]\" \"$command\"\n";
            print stderr "MEMLIMIT: $command\n";
        }
        elsif(/^TERM_OWNER:/){
            print stderr "OWNER: $command\n";
        }
    }
    close $fh_file;
}
if($isoutput==0){
    print "$count jobs checked\n";
    print stderr "All jobs finished properly!\n";
}
else{
    print stderr "\n";
}
close $fh_outfile;
unlink $outfile if($isoutput==0);
