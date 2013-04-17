#!/usr/bin/perl
use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","o=s","l:i","e:i","help");
if((!defined $opts{i})){
    &Usage();
}

my $blastfile=$opts{i};
my $outputfile=(defined $opts{o}) ? $opts{o} : "$blastfile.anno";
my $intronlength=(defined $opts{l}) ? $opts{l} : 10000;
my $extendlength=(defined $opts{e}) ? $opts{e} : 0;

open(my $fh_blastfile,$blastfile) || die "Can't open blast file\n";

my $queryID;
my $subID;
my @HSP;
my $isMergeHit=0;
my %hitinfors;

while(<$fh_blastfile>){
    chomp();
    my @lines=split(/\t/,$_);
    if($queryID ne $lines[0]){
        if($isMergeHit==0 && $queryID ne ""){
            &FormatData(@HSP);
        }
        $queryID=$lines[0];
        @HSP=();
        push(@HSP,$_);
        $subID=$lines[1];
        $isMergeHit=0;
    }
    else{
        if($subID ne $lines[1]){
            &FormatData(@HSP);
            $isMergeHit=1;
        }
        else{
            push(@HSP,$_);
        }
    }
}
&FormatData(@HSP);

open(my $fh_outfile,">$outputfile");
foreach my $key (keys %hitinfors){
    my %hitinfo=%{$hitinfors{$key}};
    my @exons=@{$hitinfo{'exons'}};
    print $fh_outfile "$hitinfo{'geneid'}\t$hitinfo{'chrom'}\t$hitinfo{'strand'}\t";
    my $chrS=$hitinfo{'chrS'}-$extendlength;
    my $chrE=$hitinfo{'chrE'}+$extendlength;
    $chrS=0 if($chrS<0);
    print $fh_outfile "$chrS\t$chrE\t$hitinfo{'exonnum'}";
    my $starts;
    my $lengths;
    foreach my $exon (@exons){
        my %exoninfo=%{$exon};
        my $exonlength=$exoninfo{'te'}-$exoninfo{'ts'};
        if($starts eq ""){
            $starts=$exoninfo{'ts'};
            $lengths=$exonlength;
        }
        else{
            $starts="$starts,$exoninfo{'ts'}";
            $lengths="$lengths,$exonlength";
        }
    }
    print $fh_outfile "\t$starts\t$lengths\n";
}
close $fh_outfile;
exit(1);

sub FormatData(){
    my @hits=@_;
    my @infos=split(/\t/,$hits[0]);
    my $hitcenter=($infos[8]+$infos[9])/2;
    my $strand=($infos[8]<$infos[9]) ? "+": "-";
    my @exons;
    foreach my $hit(@hits){
        my %exoninfo;
        my @hitsinfo=split(/\t/,$hit);
        my $str=($hitsinfo[8]<$hitsinfo[9]) ? "+": "-";
        next if($str ne $strand);
        my $exoncenter=($hitsinfo[8]+$hitsinfo[9])/2;
        next if(abs($exoncenter-$hitcenter)>$intronlength);
        $exoninfo{'qs'}=$hitsinfo[6];
        $exoninfo{'qe'}=$hitsinfo[7];
        if($strand eq "+"){
            $exoninfo{'ts'}=$hitsinfo[8];
            $exoninfo{'te'}=$hitsinfo[9];
        }
        else{
            $exoninfo{'ts'}=$hitsinfo[9];
            $exoninfo{'te'}=$hitsinfo[8];    
        }
        push(@exons,\%exoninfo);
    }
    my @sortexons=&sortexons(@exons);
    my ($chrS,$chrE,$exonnum)=&chromregion(@sortexons);
    $hitinfors{$infos[0]}->{'chrom'}=$infos[1];
    $hitinfors{$infos[0]}->{'strand'}=$strand;
    $hitinfors{$infos[0]}->{'geneid'}=$infos[0];
    $hitinfors{$infos[0]}->{'chrS'}=$chrS;
    $hitinfors{$infos[0]}->{'chrE'}=$chrE;
    $hitinfors{$infos[0]}->{'exonnum'}=$exonnum;
    $hitinfors{$infos[0]}->{'exons'}=\@sortexons;
}

sub sortexons{
    my @exons=@_;
    my @sorted;
    while(scalar(@exons)>0){
        my $minindex;
        my $minstart;
        my %exoninfo=%{$exons[0]};
        $minstart=$exoninfo{'ts'};
        for(my $i=1;$i<@exons;$i++){
            my %exoninfoN=%{$exons[$i]};
            my $qstart=$exoninfoN{'ts'};
            if($qstart<$minstart){
                $minstart=$qstart;
                $minindex=$i;
            }
        }
        push(@sorted,$exons[$minindex]);
        splice(@exons,$minindex,1);
    }
    return @sorted;
}

sub chromregion{
    my @exons=@_;
    my %exon1=%{$exons[0]};
    my $chrS=$exon1{'ts'};
    my $arraylength=scalar(@exons);
    my %exonE=%{$exons[$arraylength-1]};
    my $chrE=$exonE{'te'};
    for(my $i=0;$i<@exons;$i++){
	my $start=${$exons[$i]}{'ts'};
	my $end=${$exons[$i]}{'te'};
	if($start<$chrS){
	    $chrS=$start;
	}
	if($end>$chrE){
	    $chrE=$end;
	}
    }
    return ($chrS,$chrE,$arraylength);
}


sub Usage(){
  print << "    Usage";

	Usage:  $0 

	<options>
		-i     Input blast file (m8 format)
		-o     Output file (default: inputfile.mrg)
		-l     Maximum intron length (default: 10000)
                -e     Extend length of the gene region (default: 0)
    Usage

	exit(0);
};