#!/usr/bin/perl

#===============================================================================
#
#         FILE: SV_VCFMerge.pl
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
#      VERSION: 1.1
#      CREATED: 2013-01-08
#     REVISION: 2013-02-05
#===============================================================================

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","t=s","m:i","o:i","a:i","b:i");
if(!defined($opts{i}) || !defined($opts{t})){
    &Usage();
    exit(1);
}

my $VCFfile=$opts{i};
my $SVtype=$opts{t};
my $Outwardlen=$opts{a};
my $Inwardlen=$opts{b};
my $MergeMode=$opts{m};
my $Overlap=$opts{o};

$Outwardlen ||= 50;
$Inwardlen ||= 100;
$MergeMode ||= 0;
$Overlap ||= 0.5;

my @FileName=split(/\./,$VCFfile);
pop(@FileName);
push(@FileName,"mrg.vcf");

my $Outfile=join(".",@FileName);
my $lastchrom="";
my $lastS;
my $lastE;
my @brks;

open(my $fh_filein,"$VCFfile") || die "$!";
open(my $fh_out,">$Outfile");

while(<$fh_filein>){
    if(/^#/){
        print $fh_out "$_";
        next;
    }
    chomp();
    my @lines=split(/\t/,$_);
    my @infos=split(/;/,$lines[7]);
    my $endsite=0;
    for my $info(@infos){
        if($info=~/END=(\d+)/){
            $endsite=$1;
            last;
        }
    }
    if($lines[0] ne $lastchrom){
        if($lastchrom ne ""){
            &Output_Brks(@brks);   
        }
        @brks=();
        my %info;
        $info{'s'}=$lines[1];
        $info{'e'}=$endsite;
        push(@brks,\%info);
        $lastchrom=$lines[0];
        $lastS=$lines[1];
        $lastE=$endsite;
    }
    else{
        my $isoverlap=0;
        if($MergeMode==0){
            $isoverlap=&check_overlap($lines[1],$endsite);
        }
        else{
            $isoverlap=&check_break($lines[1],$endsite);
        }
        if($isoverlap == 0){
            my %info;
            $info{'s'}=$lines[1];
            $info{'e'}=$endsite;
            push(@brks,\%info);
        }
    }
}
&Output_Brks(@brks);
close $fh_filein;
close $fh_out;

exit(0);

sub Output_Brks{
    for(my $i=0;$i<@brks;$i++){
        my %brkinfo=%{$brks[$i]};
        my $length=$brkinfo{'e'}-$brkinfo{'s'}+1;
        if($SVtype eq "DEL"){
            $length="-$length";
        }
        print $fh_out "$lastchrom\t$brkinfo{'s'}\t.\t.\t\<$SVtype\>\t.\tPASS\tEND=$brkinfo{'e'}\;SVTYPE=$SVtype\;SVLEN=$length\n";
    }
}


sub check_overlap{
    my ($start,$end)=@_;
    for(my $i=0;$i<@brks;$i++){
        my %brkinfo=%{$brks[$i]};
        my $svlength=$brkinfo{'e'}-$brkinfo{'s'}+1;
        next if($brkinfo{'e'} < $start);
        if($brkinfo{'s'} <= $start && $brkinfo{'e'} > $start){
            if($brkinfo{'e'} >= $end){
                my $overlength=$end-$start;
                return 1 if($overlength/$svlength >= $Overlap);
            }
            else{
                my $overlength=$brkinfo{'e'}-$start;
                my $newsvlength=$end-$start+1;
                if($overlength/$svlength >= $Overlap && $overlength/$newsvlength >= $Overlap){
                    $brks[$i]->{'e'}=$end;
                    return 1;
                }
            }
        }
        elsif($brkinfo{'s'} >= $start && $brkinfo{'s'} <=$end){
            if($brkinfo{'e'} <= $end){
                my $newsvlength=$end-$start+1;
                if($svlength/$newsvlength >= $Overlap){
                    $brks[$i]->{'s'}=$start;
                    $brks[$i]->{'e'}=$end;
                    return 1;
                }
            }
            else{
                my $overlength=$end-$brkinfo{'s'};
                my $newsvlength=$end-$start+1;
                if($overlength/$svlength >= $Overlap && $overlength/$newsvlength >= $Overlap){
                    $brks[$i]->{'s'}=$start;
                    return 1;
                }
            }
        }
        last if($brkinfo{'s'} > $end);
    }
    return 0;
}

sub check_break{
	my($start,$end)=@_;
	for(my $i=0;$i<@brks;$i++){
		my %brkinfo=%{$brks[$i]};
		next if($brkinfo{'e'} < $start);
		return 1 if($start >= $brkinfo{'s'}-$Outwardlen && $start <= $brkinfo{'s'}+$Inwardlen && $end >= $brkinfo{'e'}-$Inwardlen && $end <= $brkinfo{'e'}+$Outwardlen);
		last if($brkinfo{'s'} > $end);
	}
	return 0;
}

sub Usage {#help subprogram
    print << "    Usage";

	Usage: -i <VCF_File> -t <SV_Type:INS,DEL,DUP,INV> [options]

        Options:
                 -m     Filter Mode (default: 0-overlap , 1-breakpoint)
                 
                 -o     Overlap ratio (default: 0.5)
        
                 -a     Breakpoint outward flanking length (default:50)
        
                 -b     Breakpoint inward length (default: 100)

    Usage

    exit(0);
};
