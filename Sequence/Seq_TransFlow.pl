#!/usr/bin/perl
##This program is used to analysis trancriptome informations based on reference genome
use strict;

my $version="1.0 Alvin Chen 2010-2-26";

use Getopt::Long;
my %opts;
#c configure file, i input file, s score cut for tagreads, a R space file, d Database to match,t Tag length, r redundance model
GetOptions(\%opts,"i=s","d=s","c:s","s:s","r:i","t:i","tl:i","a:s","cl:i","f:i","g:i","m:i","rt:s","l:s","td:s","st:i","help");
if((!defined $opts{i})||(!defined $opts{d})){
    &Usage();
}
##parameters
my $timestart=time;

my $parInputFile=$opts{i};
my $parDatabase=$opts{d};
my $par
my $parConffile=(defined $opts{c}) ? $opts{c} : "config.txt";
my $parScorecut=(defined $opts{s}) ? $opts{s} : 20;
my $parRdmode=(defined $opts{r}) ? $opts{r} : 1;
my $parTaglength=(defined $opts{tl}) ? $opts{tl} : 17;
my $parCountCutoff=(defined $opts{f}) ? $opts{f} : 5;
my $parCleanfile=(defined $opts{e}) ? $opts{e} : 0;
my $parSegnumber=(defined $opts{g}) ? $opts{g} : 30000;
my $parMismatch=(defined $opts{m}) ? $opts{m} : 0;
my $parTopnumber=(defined $opts{t}) ? $opts{t} : 200;
my $parRTemp=(defined $opts{rt}) ? $opts{rt} : "temple.R";
my $parTargetDir=(defined $opts{td}) ? $opts{td} : "$parSourcefile.dir";
my $parLogfile=(defined $opts{l}) ? $opts{l} : "$parSourcefile.$parTargetDir.log";
my $parClearFile=(defined $opts{cl}) ? $opts{cl} :0;
my $parSorttype=(defined $opts{st}) ? $opts{st} : 0;
