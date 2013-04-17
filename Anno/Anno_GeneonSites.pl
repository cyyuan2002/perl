#!/usr/bin/perl

#===============================================================================
#
#         FILE: 
#
#        USAGE:
#
#  DESCRIPTION:
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
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","g=s","f=i");
if(! defined($opts{i}) || ! defined($opts{g}) || ! defined($opts{f})){
    die "Usage:$0 -i <Site_file> -g <GFF_file> -f <Flanking_length>\n";
}

my $parInputFile=$opts{i};
my $parGffFile=$opts{g};
my $parFlankLen=$opts{f};
