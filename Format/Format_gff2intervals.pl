#!/usr/bin/perl

#===============================================================================
#
#         FILE: 
#
#        USAGE:
#
#  DESCRIPTION:This script is used to create interval files of intergenic region
#              and introns using GFF3 file for GATK
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

sub Sort_Exons{
    my @exons=@_;
    my @sorted;
    @sorted = sort {$a->{'start'} <=> $b->{'start'} } @exons;
    return @sorted;
}
sub Get_ChrGenes{
    my ($geneAnno,$chr)=@_;
    my @chrgenes;
    my %anno=%{$geneAnno};
    foreach my $key(keys %anno){
        if($anno{$key}->{'chrom'} eq  $chr){
            push (@chrgenes,$anno{$key});
        }
    }
    return @chrgenes;
}


sub Sort_Chroms{
    my ($geneAnno,$chrom)=@_;
    my @chroms=sort (keys %{$chrom});
    my %Chrom_Anno;
    foreach my $chr (@chroms){
        my @chrgenes = Get_ChrGenes($geneAnno,$chr);
        my @sorted = sort {$a->{'start'} <=> $b->{'start'}} @chrgenes;
        $Chrom_Anno{$chr}=\@sorted;
    }
    return \%Chrom_Anno;
}

sub Get_UTR{
    my ($cds,$mRNA_s,$mRNA_e)=@_;
    my ($UTRA,$UTRB);
    my ($cdsF,$cdsR);

    if (@{$cds} > 1) {
        $cdsF=shift(@{$cds});
        $cdsR=pop(@{$cds});
    }
    else{
        $cdsF=@{$cds}[0];
        $cdsR=@{$cds}[0];
    }

    my $Len_UTRa=$cdsF->{'start'} - $mRNA_s;
    if ($Len_UTRa > 0) {
        $UTRA = $cdsF->{'start'}-1;
    }
    else{
        $UTRA="NA";
    }

    my $Len_UTRb=$mRNA_e - $cdsR->{'end'};
    if($Len_UTRb > 0){
        $UTRB = $cdsR->{'end'}+1;
    }
    else{
        $UTRB="NA";
    }

    return ($UTRA, $UTRB);
}
sub GFF_Reader{
#----Data structure---------
#    GeneID
#       name
#	chrom
#	start
#	end
#	strand
#	transcript
#               [0]->ID
#               [0]->parent
#               [0]->name
#		[0]->start
#		[0]->end
#		[0]->5UTR
#		[0]->3UTR
#		[0]->exons
#			[0]->start
#			[0]->end
#			[1]->start
#			[1]->end
#		[1]->start
#		[1]->end
#		[1]->5UTR
#		[1]->3UTR
#		[1]->exons
#			[0]->start
#			[0]->end
#			[1]->start
#			[1]->end
#-------------------------------
    my ($GFF_file,$mode)=@_;  ## mode 1: genome mode, mode 2=gene mode
    my %geneAnno;
    my %chroms;
    open(my $gff_file,$GFF_file) or die "Can't open file: $GFF_file";
    my $lastGN;
    my $lastTN;
    my $lastStrand;
    my @transcripts;
    my %transcript;
    my @exons;
    my @cds;
    while (<$gff_file>) {
        chomp();
        next if ($_ eq "");
        my @lines=split(/\t/);
        $chroms{$lines[0]}=1;
        my ($start,$end);
        if ($lines[3] < $lines[4]) {
            $start=$lines[3];
            $end=$lines[4];
        }
        else{
            $start=$lines[4];
            $end=$lines[3];
        }
        if ($lines[2] eq 'gene') {
            if ($lastGN ne "") {
                @exons=Sort_Exons(@exons);
                @cds=Sort_Exons(@cds);
                my @tmpexons=@exons;
                my @tmpcds=@cds;
                $transcript{'cds'}=\@tmpcds;
                $transcript{'exons'}=\@tmpexons;
                my ($UTRA,$UTRB)=Get_UTR(\@cds,$transcript{'start'},$transcript{'end'});
                $transcript{'UTRA'}=$UTRA;
                $transcript{'UTRB'}=$UTRB;
                my %tmptrans=%transcript;
                push(@transcripts,\%tmptrans);
                my @tmptrans=@transcripts;
                $geneAnno{$lastGN}->{'transcripts'}=\@tmptrans;
            }
            my @geneinfo=split(/\;/,$lines[8]);
            my ($GeneID,$GeneName);
            foreach my $info(@geneinfo){
                if ($info=~/ID=(.*)/) {
                    $GeneID=$1;
                }
                elsif($info=~/Name=(.*)/){
                    $GeneName=$1;
                }
            }
            $geneAnno{$GeneID}->{'ID'}=$GeneID;
            $geneAnno{$GeneID}->{'name'}=$GeneName;
            $geneAnno{$GeneID}->{'chrom'}=$lines[0];
            $geneAnno{$GeneID}->{'start'}=$start;
            $geneAnno{$GeneID}->{'end'}=$end;
            $geneAnno{$GeneID}->{'strand'}=$lines[6];
            $lastGN=$GeneID;
            $lastStrand=$lines[6];
            @transcripts=();
            %transcript=();
            @exons=();
            @cds=();
        }
        elsif($lines[2] eq 'mRNA'){
            if (exists($transcript{'start'})) {
                @exons=Sort_Exons(@exons);
                @cds=Sort_Exons(@cds);
                my @tmpexons=@exons;
                my @tmpcds=@cds;
                $transcript{'cds'}=\@tmpcds;
                $transcript{'exons'}=\@tmpexons;
                my ($UTRA,$UTRB)=Get_UTR(\@cds,$transcript{'start'},$transcript{'end'});
                $transcript{'UTRA'}=$UTRA;
                $transcript{'UTRB'}=$UTRB;
                my %tmptrans=%transcript;
                push(@transcripts,\%tmptrans);
                @exons=();
                @cds=();
                %transcript=();
            }
            my @geneinfo=split(/\;/,$lines[8]);
            my ($GeneID,$GeneName,$ParentID);
            foreach my $info(@geneinfo){
                if ($info=~/ID=(.*)/) {
                    $GeneID=$1;
                }
                elsif($info=~/Name=(.*)/){
                    $GeneName=$1;
                }
                elsif($info=~/Parent=(.*)/){
                    $ParentID=$1;
                }
            }
            $transcript{'ID'}=$GeneID;
            $transcript{'name'}=$GeneName;
            $transcript{'parent'}=$ParentID;
            $transcript{'start'}=$start;
            $transcript{'end'}=$end;
        }
        elsif($lines[2] eq 'exon'){
            my %exon;
            $exon{'start'}=$start;
            $exon{'end'}=$end;
            push(@exons,\%exon);
        }
        elsif($lines[2] eq 'CDS'){
            my %CDS;
            $CDS{'start'}=$start;
            $CDS{'end'}=$end;
            push(@cds,\%CDS);
        }
    }
    {
        @exons=Sort_Exons(@exons);
        @cds=Sort_Exons(@cds);
        my @tmpexons=@exons;
        my @tmpcds=@cds;
        $transcript{'cds'}=\@tmpcds;
        $transcript{'exons'}=\@tmpexons;
        my ($UTRA,$UTRB)=Get_UTR(\@cds,$transcript{'start'},$transcript{'end'});
        $transcript{'UTRA'}=$UTRA;
        $transcript{'UTRB'}=$UTRB;
        my %tmptrans=%transcript;
        push(@transcripts,\%tmptrans);
        my @tmptrans=@transcripts;
        $geneAnno{$lastGN}->{'transcripts'}=\@tmptrans;
    }
    my $chr_anno=Sort_Chroms(\%geneAnno,\%chroms);
    if ($mode ==1) {
        return $chr_anno;
    }
    elsif($mode ==2){
        return \%geneAnno;
    }
}

#----------------Main start here-------

my ($GffFile,$ChromLenFile)=@ARGV;
if (@ARGV<2) {
    die "Usage:$0 <GFF_File> <Chrom_Length_File>\n";
}
my %chrom_len;

open(my $fh_chrom,"<",$ChromLenFile) or die "Can't open file: $ChromLenFile";
while (<$fh_chrom>) {
    chomp();
    my @lines=split(/\t/,$_);
    $chrom_len{$lines[0]}=$lines[1];
}
close $fh_chrom;

if (!-e($GffFile)) {
    die "Can't open file: $GffFile\n";
}
my $refGffAnno=GFF_Reader($GffFile,1);

foreach my $key(sort keys %chrom_len){
    my @chrominfo=@{$refGffAnno->{$key}};
    my $startpos=1;
    for(my $i=0;$i<@chrominfo;$i++){
        my $refexons=$chrominfo[$i]->{'transcripts'}->[0]->{'exons'};
        my @exons=@{$refexons};
        for(my $j=0;$j<@exons;$j++){
            my $endpos=$exons[$j]->{'start'}-1;
            print "$key:$startpos-$endpos\n";
            $startpos=$exons[$j]->{'end'}+1;
        }
    }
    print "$key:$startpos-$chrom_len{$key}\n";
}

exit(0);