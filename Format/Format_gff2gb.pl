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
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;

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
                my @tmpexons=@exons;
                my @tmpcds=@cds;
                @exons=Sort_Exons(@exons);
                @cds=Sort_Exons(@cds);
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
                my @tmpexons=@exons;
                my @tmpcds=@cds;
                @exons=Sort_Exons(@exons);
                @cds=Sort_Exons(@cds);
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
        my @tmpexons=@exons;
        my @tmpcds=@cds;
        @exons=Sort_Exons(@exons);
        @cds=Sort_Exons(@cds);
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

sub gffinfo2gb{
    my($gID,$refGffAnno,$RefGenome,$Flanking_len)=@_;
    my %GffAnno=%{$refGffAnno};
    if (exists($GffAnno{$gID})) {
            my $Chrom=$GffAnno{$gID}->{'chrom'};
            my $gS=$GffAnno{$gID}->{'start'};
            my $gE=$GffAnno{$gID}->{'end'};
            my $rS=$gS-$Flanking_len;
            my $rE=$gE+$Flanking_len;
            my $strand=$GffAnno{$gID}->{'strand'};
            my $seqregion="$Chrom:$rS-$rE";
            my $seq=$RefGenome->seq($seqregion);
            if ($strand eq "-") {
                $seq=~tr/atgcATGC/tacgTACG/;
                $seq=reverse($seq);
            }
            
            my $seq_obj = Bio::Seq->new(-seq => $seq, display_id => "$gID");
            my ($refgS,$refgE);
            if ($strand eq "+") {
                $refgS=$gS-$rS+1;
                $refgE=$gE-$rS+1;
            }
            else{
                $refgS=$rE-$gE+1;
                $refgE=$rE-$gS+1;
            }
            my $refLen=$rE-$rS+1;
            my $seqfeat1=new Bio::SeqFeature::Generic(-start => 1,
                                                      -end => $refLen,
                                                      );
            $seq_obj->add_SeqFeature($seqfeat1);
            my $seqfeat2=new Bio::SeqFeature::Generic(-start => $refgS,
                                                     -end => $refgE,
                                                     -seq_id => $gID,
                                                     -primary_tag => 'gene');
            $seq_obj->add_SeqFeature($seqfeat2);
            my @transcripts=@{$GffAnno{$gID}->{'transcripts'}};            
            foreach my $transinfo(@transcripts){
                my ($tranS,$tranE);
                if ($strand eq "+") {
                    $tranS=$transinfo->{'start'}-$rS+1;
                    $tranE=$transinfo->{'end'}-$rS+1;
                }
                else{
                    $tranS=$rE-$transinfo->{'end'}+1;
                    $tranE=$rE-$transinfo->{'start'}+1;
                }
                my @exons=@{$transinfo->{'exons'}};
                my $regioninfo;
                foreach my $exon(@exons){
                    my($refS,$refE);
                    if ($strand eq "+") {
                        $refS=$exon->{'start'}-$rS+1;
                        $refE=$exon->{'end'}-$rS+1;
                    }
                    else{
                        $refS=$rE-$exon->{'end'}+1;
                        $refE=$rE-$exon->{'start'}+1;
                    }
                    if ($regioninfo eq "") {
                        $regioninfo="$refS..$refE";
                    }
                    else{
                        $regioninfo.=",$refS..$refE";
                    }
                }
                my $mRNAfeat= new Bio::SeqFeature::Generic(-start => $tranS,
                                                          -end => $tranE,
                                                          -phrase => $regioninfo,
                                                          -primary_tag => 'mRNA'
                                                          );
                
                $seq_obj->add_SeqFeature($mRNAfeat);
                my @cdss=@{$transinfo->{'cds'}};
                $regioninfo="";
                my $cdsS=0;
                my $cdsE=0;
                foreach my $cds(@cdss){
                    my($refS,$refE);
                    if ($strand eq "+") {
                        $refS=$cds->{'start'}-$rS+1;
                        $refE=$cds->{'end'}-$rS+1;
                    }
                    else{
                        $refS=$rE-$cds->{'end'}+1;
                        $refE=$rE-$cds->{'start'}+1;
                    }
                    if ($cdsS == 0) {
                        $cdsS=$refS;
                    }
                    else{
                        if ($refS < $cdsS) {
                            $cdsS=$refS;
                        }
                    }
                    if ($refE > $cdsE) {
                        $cdsE = $refE;
                    }
                    
                    
                    if ($regioninfo eq "") {
                        $regioninfo="$refS..$refE";
                    }
                    else{
                        $regioninfo.=",$refS..$refE";
                    }
                }
                my $cdsfeat= new Bio::SeqFeature::Generic(-start => $cdsS,
                                      -end => $cdsE,
                                      -phrase => $regioninfo,
                                      -primary_tag => 'cds'
                                      );
                $seq_obj->add_SeqFeature($cdsfeat);
                return ($seq_obj);
            }            
        }
        else{
            print stderr "Opps: Can't find $gID\n";
            return "";
        }
}

my %opts;
my $usage="$0 <-s sequence_file> <-g gff_file> [options]\n";
$usage.="       -f flanking_length, default 0\n";
$usage.="       -l file of output gene list\n";

GetOptions(\%opts,'s=s','g=s','f:i','l:s');
if (!exists($opts{s}) || !exists($opts{g})) {
    die "$usage";    
}

my $SeqFile=$opts{'s'};
my $GffFile=$opts{'g'};
my $Flanking_len=$opts{'f'};
my $GeneListFile=$opts{'l'};

$Flanking_len||=0;

if (!-e($SeqFile)) {
    die "Can't open file: $SeqFile\n";
}
if (!-e($GffFile)) {
    die "Can't open file: $GffFile\n";
}

my $RefGenome=Bio::DB::Fasta->new($SeqFile);
my $refGffAnno=GFF_Reader($GffFile,2);

if ($GeneListFile ne "") {
    open(my $fh_glfile,"$GeneListFile") or die "Can't open file: $GeneListFile";
    my $outfile="$GeneListFile.gb";
    my $io = Bio::SeqIO->new(-format => "genbank", -file=>">$outfile");
    while (<$fh_glfile>) {
        chomp();
        my $gID=$_;
        my $gbseq = gffinfo2gb($gID,$refGffAnno,$RefGenome,$Flanking_len);
        if ($gbseq ne "") {
            $io->write_seq($gbseq);
        }
    }
    close $fh_glfile;
}
else{
    my $outfile="$GffFile.gb";
    my $io = Bio::SeqIO->new(-format => "genbank", -file=>">$outfile");
    foreach my $gID(keys %{$refGffAnno}){
        my $gbseq=gffinfo2gb($gID,$refGffAnno,$RefGenome,$Flanking_len);
        if ($gbseq ne "") {
            $io->write_seq($gbseq);
        }
    }
}

exit(0);

