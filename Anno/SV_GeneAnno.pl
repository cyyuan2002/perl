#!/usr/bin/perl
#===============================================================================
#
#         FILE:SV_GeneAnno.pl
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
#      VERSION:
#      CREATED:
#     REVISION:
#===============================================================================


#-------Breakpoint File Format---------
#supercont2.1	339475	supercont2.1	350400	DEL
#supercont2.1	485505	supercont2.1	485583	DEL
#supercont2.1	485505	supercont2.1	485583	INV
#-------------------------------------

use strict;

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

sub Sort_Exons{
    my @exons=@_;
    my @sorted;
    @sorted = sort {$a->{'start'} <=> $b->{'start'} } @exons;
    return @sorted;
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
    my $GFF_file=shift;
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

    return $chr_anno;
}

sub Search_GeneRegion{
    my ($BPSite,$chrom_info,$geneindex)=@_;
    my @Chrom_Anno=@{$chrom_info};
    my $geneinfo=$Chrom_Anno[$geneindex];
    my $gene=$geneinfo->{'ID'};
    my $genestrand=$geneinfo->{'strand'};
    my $transcript=$geneinfo->{'transcripts'}->[0];
    my @exons=@{$transcript->{'exons'}};
    my $lastend;
    my $outtext;
    for(my $i=0;$i<@exons;$i++){
        if ($i > 0) {
            if ($lastend <= $BPSite && $exons[$i]->{'start'} > $BPSite) {
                my $intro_count=$i;
                if ($geneinfo->{'strand'} eq "-") {
                    $intro_count=scalar(@exons)-$intro_count;
                }
                $outtext = "intro_$intro_count -- $gene($genestrand)\t";
                last;
            }
        }
        if ($exons[$i]->{'start'} <= $BPSite && $exons[$i]->{'end'} >=$BPSite) {
            if ($BPSite <= $transcript->{'UTRA'}) {
                if ($geneinfo->{'strand'} eq "+") {
                    $outtext = "5UTR -- $gene($genestrand)\t";
                }
                else{
                    $outtext = "3UTR -- $gene($genestrand)\t";
                }
            }
            elsif ($BPSite >= $transcript->{'UTRB'}) {
                if ($geneinfo->{'strand'} eq "+") {
                    $outtext = "3UTR -- $gene($genestrand)\t";
                }
                else{
                    $outtext = "5UTR -- $gene($genestrand)\t";
                }
            }
            else{
                my $exon_count;
                if($geneinfo->{'strand'} eq "-"){
                    $exon_count=scalar(@exons)-$i;
                }
                else{
                    $exon_count=$i+1;
                }
                $outtext = "exon_$exon_count -- $gene($genestrand)\t";
            }
            last;
        }
        $lastend=$exons[$i]->{'end'};
    }

    return $outtext;
}

sub Search_MidRegion{
    my ($BPleft,$BPright,$chrom_info)=@_;
    my @Chrom_Anno=@{$chrom_info};
    my $genes;
    for(my $i=0;$i<@Chrom_Anno;$i++){
        my $geneinfo=$Chrom_Anno[$i];
        if ($geneinfo->{'start'} >= $BPright) {
            last;
        }
        elsif(($geneinfo->{'start'} <= $BPleft && $geneinfo->{'end'} >= $BPleft) ||
              ($geneinfo->{'start'} <= $BPright && $geneinfo->{'end'} >=$BPright) ||
              ($BPleft <= $geneinfo->{'start'} && $BPright >= $geneinfo->{'end'})){
            if ($genes ne "") {
                $genes=$genes.",".$geneinfo->{'ID'};
            }
            else{
                $genes=$geneinfo->{'ID'};
            }
        }
    }
    $genes="NA" if ($genes eq "");
    return $genes;
}

sub Search_BPLocus{
    my ($bpsite,$chrom_info,$isleft)=@_;
    my @Chrom_Anno=@{$chrom_info};
    my $lastEnd=0;
    my $isBPfound=0;
    my $typeBP="";
    my $indexBP="";
    for(my $i=0;$i<@Chrom_Anno;$i++){
        my $geneinfo=$Chrom_Anno[$i];
        if ($isBPfound==0) {
            if ($lastEnd <= $bpsite && $geneinfo->{'start'} > $bpsite) {
                $isBPfound=1;
                $typeBP='I'; ## "intergenic"
                if ($isleft == 1) {
                    $indexBP=$i-1;
                }
                else{
                    $indexBP=$i;
                }
                last;
            }
            elsif ($geneinfo->{'start'} <= $bpsite && $geneinfo->{'end'} >= $bpsite) {
                $isBPfound=1;
                $typeBP='G'; ## "generegion"
                $indexBP=$i;
                last;
            }
        }
        $lastEnd=$geneinfo->{'end'};
    }
    return ($isBPfound,$typeBP,$indexBP);
}

#-----------------------Main------------------------

my ($parGFFFile,$parBpFile)=@ARGV;
if (@ARGV<2) {
    die "Usage: <GFF3_File> <BreakPoint_File>\n";
}

my $Gff_info=GFF_Reader($parGFFFile);

open(my $fh_BpFile,"$parBpFile") or die "Can't open file: $parBpFile";
while (<$fh_BpFile>) {
    chomp();
    my @lines=split(/\t/,$_);
    my @Chrom_Anno=@{$Gff_info->{$lines[0]}};
    my $bpleft=$lines[1];
    my $bpright=$lines[3];
    print "$_\t";
    if ($lines[4] eq 'DEL' || $lines[4] eq 'DUP') {
        my($foundBPLeft,$typeBPLeft,$indexBPLeft)=Search_BPLocus($bpleft,\@Chrom_Anno,1);
        my($foundBPRight,$typeBPRight,$indexBPRight)=Search_BPLocus($bpright,\@Chrom_Anno,0);

        if ($foundBPLeft==1) {
            if ($typeBPLeft eq "I") {
                if ($indexBPLeft >=0) {
                    my $leftlength=$bpleft-$Chrom_Anno[$indexBPLeft]->{'end'};
                    my $leftgene=$Chrom_Anno[$indexBPLeft]->{'ID'};
                    my $leftgenestrand=$Chrom_Anno[$indexBPLeft]->{'strand'};
                    print "LEFT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand)\t";
                }
                else{
                    print "LEFT:intergenic -- pre_gene NA\t";
                }
            }
            else{
                my $outText=Search_GeneRegion($bpleft,\@Chrom_Anno,$indexBPLeft);
                print "LEFT:$outText\t";
            }
            if ($foundBPRight==1) {
                if ($typeBPRight eq "I") {
                    my $rightlength=$Chrom_Anno[$indexBPRight]->{'start'}-$bpright;
                    my $rightgene=$Chrom_Anno[$indexBPRight]->{'ID'};
                    my $rightgenestrand=$Chrom_Anno[$indexBPRight]->{'strand'};
                    print "RIGHT:intergenic -- next_gene($rightlength) $rightgene($rightgenestrand)\t";
                }
                else{
                    my $outText=Search_GeneRegion($bpright,\@Chrom_Anno,$indexBPRight);
                    print "RIGHT:$outText\t";
                }
            }
            else{
                print "RIGHT:intergenic -- next_gene NA\t";
            }

            ##MID region of the SV
            my $OutMid=Search_MidRegion($bpleft,$bpright,\@Chrom_Anno);
            print "MID:$OutMid\n";
        }
        else{ ## 3' of the chromsome
            my $leftgeneinfo=$Chrom_Anno[scalar(@Chrom_Anno)-1];
            my $leftlength=$bpleft-$leftgeneinfo->{'end'};
            my $leftgene=$leftgeneinfo->{'ID'};
            my $leftgenestrand=$leftgeneinfo->{'strand'};
            print "LEFT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand)\t";
            print "RIGHT:intergenic -- next_gene NA\t";
            print "MID:NA\n";
        }
    }
    elsif($lines[4] eq 'INV' || $lines[4] eq 'CTX'){
        my($foundBPLeft,$typeBPLeft,$indexBPLeft)=Search_BPLocus($bpleft,\@Chrom_Anno,1);
        if ($foundBPLeft == 1) {
            if ($typeBPLeft eq "I") {
                if ($indexBPLeft >=0) {
                    my $leftlength=$bpleft-$Chrom_Anno[$indexBPLeft]->{'end'};
                    my $leftgene=$Chrom_Anno[$indexBPLeft]->{'ID'};
                    my $leftgenestrand=$Chrom_Anno[$indexBPLeft]->{'strand'};
                    print "LEFT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand) <|> ";
                }
                else{
                    print "LEFT:intergenic -- pre_gene NA <|> ";
                }
                my $rightlength=$Chrom_Anno[$indexBPLeft+1]->{'start'}-$bpleft;
                my $rightgene=$Chrom_Anno[$indexBPLeft+1]->{'ID'};
                my $rightgenestrand=$Chrom_Anno[$indexBPLeft+1]->{'strand'};
                print "next_gene($rightlength) $rightgene($rightgenestrand)\t";
            }
            else{
                my $outText=Search_GeneRegion($bpleft,\@Chrom_Anno,$indexBPLeft);
                print "LEFT:$outText\t";
            }
        }
        else{
            my $leftgeneinfo=$Chrom_Anno[scalar(@Chrom_Anno)-1];
            my $leftlength=$bpleft-$leftgeneinfo->{'end'};
            my $leftgene=$leftgeneinfo->{'ID'};
            my $leftgenestrand=$leftgeneinfo->{'strand'};
            print "LEFT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand) <|> next_gene NA\t";
        }

        my @Chrom_Anno2;
        if ($lines[4] eq 'CTX') {
            @Chrom_Anno2=@{$Gff_info->{$lines[2]}};
        }
        else{
            @Chrom_Anno2=@Chrom_Anno;
        }
        my ($foundBPRight,$typeBPRight,$indexBPRight)=Search_BPLocus($bpright,\@Chrom_Anno2,1);
        if ($foundBPRight == 1) {
            if ($typeBPRight eq "I") {
                if ($indexBPRight >=0) {
                    my $leftlength=$bpright-$Chrom_Anno2[$indexBPRight]->{'end'};
                    my $leftgene=$Chrom_Anno2[$indexBPRight]->{'ID'};
                    my $leftgenestrand=$Chrom_Anno2[$indexBPRight]->{'strand'};
                    print "RIGHT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand) <|> ";
                }
                else{
                    print "RIGHT:intergenic -- pre_gene NA <|> ";
                }
                my $rightlength=$Chrom_Anno2[$indexBPRight+1]->{'start'}-$bpright;
                my $rightgene=$Chrom_Anno2[$indexBPRight+1]->{'ID'};
                my $rightgenestrand=$Chrom_Anno2[$indexBPRight+1]->{'strand'};
                print "next_gene($rightlength) $rightgene($rightgenestrand)\n";
            }
            else{
                my $outText=Search_GeneRegion($bpright,\@Chrom_Anno2,$indexBPRight);
                print "RIGHT:$outText\n";
            }
        }
        else{
            my $leftgeneinfo=$Chrom_Anno2[scalar(@Chrom_Anno2)-1];
            my $leftlength=$bpright-$leftgeneinfo->{'end'};
            my $leftgene=$leftgeneinfo->{'ID'};
            my $leftgenestrand=$leftgeneinfo->{'strand'};
            print "RIGHT:intergenic -- pre_gene($leftlength) $leftgene($leftgenestrand) <|> next_gene NA\n";
        }
    }
}

