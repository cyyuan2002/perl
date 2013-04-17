#!/usr/bin/perl

#===============================================================================
#
#         FILE: SV_Mummer.pl
#
#        USAGE:
#
#  DESCRIPTION: This program is used to search inter-chromosomal translocation (CTX),
#               intra-chromosomal translcation (ITX) and inversion (INV) based on
#               the alignment result of MUMmer
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 10/31/2012
#     REVISION:
#===============================================================================

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp;
use File::Basename;

my $MUMMER="nucmer";
my $DELFLT="delta-filter";
my $SHOWCOORDS="show-coords";
my $CROSSMATCH="cross_match";

my %opts;

GetOptions(\%opts,"q:s","r:s","d:s","c:s","i:i","l:i","o:s","p:i","t:i","f:i","v:i");
if((!defined($opts{q}) && !defined($opts{r})) && (!defined($opts{d}))){
    &Usage();
    exit(1);
}

if(!defined($opts{o})){
    &Usage();
    exit(1);
}

my $parQueryFile=$opts{q};
my $parRefFile=$opts{r};
my $parOutfile=$opts{o};
my $parDeltaFile=$opts{d};
my $parCentFile=$opts{c};
my $parMinIdentity=$opts{i};
my $parMinLength=$opts{l};
my $parOverpercent=$opts{p};
my $parITXLength=$opts{t};
my $parFlankLength=$opts{f};
my $parIsVerify=$opts{v};
my $parCentroFilter=0;
my $parCrossOverMinLength=200;
my $parCrossOverMaxSub=10;
my $parMaxBreakpointsOverlap=10;

$parMinIdentity ||= 95;
$parMinLength ||= 500;
$parOverpercent ||= 0.8;
$parITXLength ||= 10000;
$parFlankLength ||= 1000;
$parIsVerify ||= 0;

my @Tempfiles;

my $parRepeatFile="$parOutfile.rep";
my $parSVFile="$parOutfile.sv";
my $parLogFile="$parOutfile.log";


if($parCentFile ne ""){
    $parCentroFilter = 1;
}

if($parIsVerify == 1){
    if($parQueryFile eq "" || $parRefFile eq ""){
        &Usage();
        exit(1);
    }
}

my %RepeatRegion;
my %RefCoords;
my %CentroRegion;
my %QuerySeqs;
my @arrayBreakpoints;
my %hashBreakpoints;

if($parDeltaFile eq ""){
    die "Can't open file $parQueryFile\n" if (!-e($parQueryFile));
    die "Can't open file $parRefFile\n" if (!-e($parRefFile));
    my $Querybase=basename($parQueryFile);
    `$MUMMER -p $Querybase $parRefFile $parQueryFile`;
    $parDeltaFile="$Querybase.delta";
}


my $CoordsFileA="$parDeltaFile.coords";
my $CoordsFileR="$parDeltaFile.flt.1.coords";
my $CoordsFileQ="$parDeltaFile.flt.2.coords";
`$SHOWCOORDS -THrl $parDeltaFile > $CoordsFileA`;
`$DELFLT -r $parDeltaFile > $parDeltaFile.flt1`;
`$DELFLT -q $parDeltaFile > $parDeltaFile.flt2`;
`$SHOWCOORDS -THrl $parDeltaFile.flt1 > $CoordsFileR`;
`$SHOWCOORDS -THql $parDeltaFile.flt2 > $CoordsFileQ`;

push(@Tempfiles,$CoordsFileA);
push(@Tempfiles,$CoordsFileR);
push(@Tempfiles,$CoordsFileQ);
push(@Tempfiles,"$parDeltaFile.flt1");
push(@Tempfiles,"$parDeltaFile.flt2");

##----Read QuerySequence-----
open (my $fh_queryseq,$parQueryFile);
my $lastseqN="";
my $seq="";
while(<$fh_queryseq>){
    chomp();
    if(/^>(\S+)/){
        if($seq ne ""){
            $QuerySeqs{$lastseqN}=$seq;
        }
        $lastseqN=$1;
        $seq="";
    }
    else{
        $seq.=$_;
    }
}
$QuerySeqs{$lastseqN}=$seq;
close $fh_queryseq;

##----Read Centromere-------
##supercont2.1    970001  1004000
##supercont2.2    866001  893000
##supercont2.3    1370501 1410000

$parCentroFilter = 1 if($parCentFile ne "");
if($parCentroFilter == 1){
    open (my $fh_centrofile,"$parCentFile") || die "Can't open file $parCentFile\n";
    while(<$fh_centrofile>){
        my @lines=split(/\t/,$_);
        $CentroRegion{$lines[0]}->{'s'}=$lines[1];
        $CentroRegion{$lines[0]}->{'e'}=$lines[2];
    }
    close $fh_centrofile;
}

##---Search for the Repeat regions in Reference-----
my $lastChrom="";
my $last_s;
my $last_e;
my $last_ml;
my $last_RS;
my $last_RE;
my $last_RC;
my @repeats;
open (my $fh_Coords, $CoordsFileA) || die "Can't open file $CoordsFileA\n";
while(<$fh_Coords>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastChrom ne $lines[9]){
        if($last_RS > 0){
            my %RP_region;
            $RP_region{'rs'}=$last_RS;
            $RP_region{'re'}=$last_RE;
            $RP_region{'rc'}=$last_RC;
            push(@repeats,\%RP_region);
        }
        if($lastChrom ne ""){
            my @tempArray=&MergeSite(@repeats);
            $RepeatRegion{$lastChrom}=\@tempArray;
        }
        $lastChrom=$lines[9];
        $last_s=$lines[0];
        $last_e=$lines[1];
        $last_ml=$lines[4];
        $last_RS=0;
        $last_RE=0;
        $last_RC=0;
        @repeats=();
        next;
    }
    if($last_RS == 0){
        if($lines[0] >= $last_s && $lines[0] < $last_e){
            my $overlength;
            if($lines[1] >= $last_e){
                $overlength=$last_e-$lines[0]+1;
            }
            else{
                $overlength=$lines[4];
            }
            if($overlength/$last_ml >= $parOverpercent && $overlength/$lines[4] >= $parOverpercent){
                $last_RS=$last_s;
                $last_RE=$lines[1];
                $last_RC=2;
            }
        }
    }
    else{
        if($lines[0] >= $last_RE){
            my %RP_region;
            $RP_region{'rs'}=$last_RS;
            $RP_region{'re'}=$last_RE;
            $RP_region{'rc'}=$last_RC;
            push(@repeats,\%RP_region);
            $last_RS=0;
            $last_RE=0;
            $last_RC=0;
        }
        else{
            if($lines[0] < $last_RE){
                if($lines[1] > $last_RE){
                    my $overlength=$last_RE-$lines[0]+1;
                    if($overlength/$lines[4] >= $parOverpercent){
                        $last_RE=$lines[1];
                        $last_RC++;
                    }
                    else{
                        my %RP_region;
                        $RP_region{'rs'}=$last_RS;
                        $RP_region{'re'}=$last_RE;
                        $RP_region{'rc'}=$last_RC;
                        push(@repeats,\%RP_region);
                        $last_RS=0;
                        $last_RE=0;
                        $last_RC=0;
                    }
                }
                else{
                    $last_RC++;
                }
            }
        }
    }
    $last_s=$lines[0];
    $last_e=$lines[1];
    $last_ml=$lines[4];
}

if($last_RS > 0){
    my %RP_region;
    $RP_region{'rs'}=$last_RS;
    $RP_region{'re'}=$last_RE;
    $RP_region{'rc'}=$last_RC;
    push(@repeats,\%RP_region);
}
@repeats=&MergeSite(@repeats);
$RepeatRegion{$lastChrom}=\@repeats;
close $fh_Coords;

open (my $fh_repeatfile,">$parRepeatFile");
foreach my $chrom(sort (keys %RepeatRegion)){
    my @reps=@{$RepeatRegion{$chrom}};
    for(my $i=0;$i<@reps;$i++){
        my %repinfo=%{$reps[$i]};
        my $replen=$repinfo{'re'}-$repinfo{'rs'}+1;
        print $fh_repeatfile "$chrom\t$repinfo{'rs'}\t$repinfo{'re'}\t$replen\t$repinfo{'rc'}\n";
    }
}
close $fh_repeatfile;
#-----------------

open (my $fh_RefCoords, $CoordsFileR) || die "Can't open file $CoordsFileR\n";
$lastChrom="";
my @hits;

while(<$fh_RefCoords>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastChrom ne $lines[9]){
        if($lastChrom ne ""){
            my @tempArray=@hits;
            $RefCoords{$lastChrom}=\@tempArray;
        }
        $lastChrom=$lines[9];
        @hits=();
    }
    my %mapinfo;
    $mapinfo{'ss'}=$lines[0];
    $mapinfo{'se'}=$lines[1];
    if($lines[2]<$lines[3]){
        $mapinfo{'qs'}=$lines[2];
        $mapinfo{'qe'}=$lines[3];
        $mapinfo{'dr'}="+";
    }
    else{
        $mapinfo{'qs'}=$lines[3];
        $mapinfo{'qe'}=$lines[2];
        $mapinfo{'dr'}="-";
    }
    $mapinfo{'rn'}=$lines[9];
    $mapinfo{'ml'}=$lines[5];
    $mapinfo{'qn'}=$lines[10];
    next if ($lines[5] < $parMinLength);
    next if ($lines[6] < $parMinIdentity);
    push(@hits,\%mapinfo);
}
$RefCoords{$lastChrom}=\@hits;
close $fh_RefCoords;

open (my $fh_logfile,">$parLogFile");

@hits=();
my $lastContig="";
open(my $fh_QueryCoords,$CoordsFileQ);
while(<$fh_QueryCoords>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastContig ne $lines[10]){
        if($lastContig ne ""){
            if(@hits>1){
                &checkhits($lastContig,\@hits);
            }
            else{
                #print "print $_\n";
            }
        }
        @hits=();
        $lastContig=$lines[10];
    }
    my %mapinfo;
    $mapinfo{'ss'}=$lines[0];
    $mapinfo{'se'}=$lines[1];
    if($lines[2]<$lines[3]){
        $mapinfo{'qs'}=$lines[2];
        $mapinfo{'qe'}=$lines[3];
        $mapinfo{'dr'}="+";
    }
    else{
        $mapinfo{'qs'}=$lines[3];
        $mapinfo{'qe'}=$lines[2];
        $mapinfo{'dr'}="-";
    }
    $mapinfo{'ml'}=$lines[5];
    $mapinfo{'rn'}=$lines[9];
    $mapinfo{'qn'}=$lines[10];
    $mapinfo{'id'}=$lines[6];
    push(@hits,\%mapinfo);
}
if(@hits>1){
    &checkhits($lastContig,\@hits);
}


open(my $fh_outfile,">$parSVFile");
foreach my $Chr (sort keys %hashBreakpoints){
    my %chrinfos=%{$hashBreakpoints{$Chr}};
    foreach my $site(sort {$a<=>$b} keys %chrinfos){
        my @points=@{$chrinfos{$site}};
        foreach my $point(@points){
            my %breakinfo=%{$point};
            print $fh_outfile "$Chr\t$site\t$breakinfo{'chr'}\t$breakinfo{'site'}\t$breakinfo{'type'}\n";
        }
    }
}

unlink @Tempfiles;
exit(0);

sub checkhits{
    my ($contigID,$hitref)=@_;
    my @hitinfors=@{$hitref};
    my $last_rn;
    my $last_ss;
    my $last_se;
    my $last_qs;
    my $last_qe;
    my $last_ml;
    my $last_id;
    my @sortedresult;
    for(my $i=0;$i<@hitinfors;$i++){ ##remove repeated mapped fregments
        my %mapinfo=%{$hitinfors[$i]};
        next if($mapinfo{'id'} < $parMinIdentity);
        next if($mapinfo{'ml'} < $parMinLength);
        if($last_rn eq ""){
            $last_rn=$mapinfo{'rn'};
            $last_ss=$mapinfo{'ss'};
            $last_se=$mapinfo{'se'};
            $last_qs=$mapinfo{'qs'};
            $last_qe=$mapinfo{'qe'};
            $last_ml=$mapinfo{'ml'};
            $last_id=$mapinfo{'id'};
            push(@sortedresult,\%mapinfo);
            next;
        }
        if($mapinfo{'qs'} < $last_qe){ ##hits overlapped
            if($mapinfo{'qe'} <= $last_qe){ ## new hits covered by last hits
                next;
            }
            else{
                my $overlength=$last_qe-$mapinfo{'qs'}+1;
                if($overlength/$last_ml >= $parOverpercent || $overlength/$mapinfo{'ml'} >= $parOverpercent){
                    my $last_mpscore=$last_id*$last_ml/100;
                    my $mpscore=$mapinfo{'id'}*$mapinfo{'ml'}/100;
                    if($last_mpscore >= $mpscore){
                        next;
                    }
                    else{
                        pop(@sortedresult);
                    }
                }
            }
        }
        $last_rn=$mapinfo{'rn'};
        $last_ss=$mapinfo{'ss'};
        $last_se=$mapinfo{'se'};
        $last_qs=$mapinfo{'qs'};
        $last_qe=$mapinfo{'qe'};
        $last_ml=$mapinfo{'ml'};
        $last_id=$mapinfo{'id'};
        push(@sortedresult,\%mapinfo);
    }

    # search info from ref_sort to decide whether to keep it or not
    my @filteredresult;
    for(my $i=0;$i<@sortedresult;$i++){
        my $iskeep=1;
        my $issearched=&SearchContig($sortedresult[$i]);
        next if($issearched==0);
        my %mapinfo=%{$sortedresult[$i]};
        my $map_rn=$mapinfo{'rn'};
        my $map_ss=$mapinfo{'ss'};
        my $map_se=$mapinfo{'se'};
        my $map_qn=$mapinfo{'qn'};
        my $map_ml=$mapinfo{'ml'};
        my $map_id=$mapinfo{'id'};
        my @RefMapinfos=@{$RefCoords{$map_rn}};
        for(my $j=0;$j<@RefMapinfos;$j++){
            my %refinfo=%{$RefMapinfos[$j]};
            my $ref_ss=$refinfo{'ss'};
            my $ref_se=$refinfo{'se'};
            my $ref_qn=$refinfo{'qn'};
            my $ref_ml=$refinfo{'ml'};
            my $ref_id=$refinfo{'id'};
            last if ($ref_ss > $map_se);
            next if($ref_qn eq $map_qn); ##ignore it self;
            if($map_ss <= $ref_se && $map_se >=$ref_se){ ##Only keep this hit when it can be searched from the ref_sort and the region is larger or equal to the ref_region
                if($map_ss >= $ref_ss){ #overlap
                    my $overlength=$ref_se-$map_ss;
                    if($overlength/$map_ml >= $parOverpercent || $overlength/$ref_ml >= $parOverpercent){
                        my $map_score=$map_ml*$map_id/100;
                        my $ref_score=$ref_ml*$ref_id/100;
                        if($ref_score > $map_score){
                            $iskeep=0;
                            last;
                        }
                    }
                }
            }
            elsif($map_ss >= $ref_ss && $map_se <= $ref_se){ ##map region is in a larger mapped region
                $iskeep=0;
                last;
            }
            elsif($map_ss <= $ref_ss && $map_se >= $ref_ss){
                if($map_se <= $ref_se){
                    my $overlength=$ref_se-$map_ss;
                    if($overlength/$map_ml >= $parOverpercent || $overlength/$ref_ml >= $parOverpercent){
                        my $map_score=$map_ml*$map_id/100;
                        my $ref_score=$ref_ml*$ref_id/100;
                        if($ref_score > $map_score){
                            $iskeep=0;
                            last;
                        }
                    }
                }
            }
        }
        next if($iskeep==0);
        push(@filteredresult, $sortedresult[$i]);
    }
    return if(@filteredresult<2);

    ##------search for the SV---------------
    ##use these variables again

    my %last_map;
    my @maphits;
    for(my $i=0;$i<@filteredresult;$i++){
        my $isRep=&CheckRepeat($filteredresult[$i]);
        next if($isRep == 1);
        if($parCentroFilter == 1){
            my $isCentro=&CheckCentro($filteredresult[$i]);
            next if($isCentro == 1);
        }
        my %mapinfo=%{$filteredresult[$i]};
        push(@maphits,$filteredresult[$i]);
        #print "$mapinfo{'rn'}\t$mapinfo{'ss'}\t$mapinfo{'se'}\t$mapinfo{'dr'}\t$mapinfo{'qs'}\t$mapinfo{'qe'}\t$contigID\n";
        if(!%last_map){ ##assign first value
            %last_map=%mapinfo;
            next;
        }
        if($last_map{'rn'} ne $mapinfo{'rn'}){##CTX
            pop(@maphits);
            if(@maphits >= 2){
                &Check_Order(\@maphits);
            }
            @maphits=();

            my $BP_ChrA;
            my $BP_SiteA;
            my $BP_ChrB;
            my $BP_SiteB;
            $BP_ChrA=$last_map{'rn'};
            $BP_ChrB=$mapinfo{'rn'};
            if($last_map{'dr'} eq "+"){
                $BP_SiteA = $last_map{'se'};
            }
            else{
                $BP_SiteA = $last_map{'ss'};
            }
            if($mapinfo{'dr'} eq "+"){
                $BP_SiteB = $mapinfo{'ss'};
            }
            else{
                $BP_SiteB = $mapinfo{'se'};
            }

            if($parIsVerify == 1){
                my $contig_s=$last_map{'qe'}-$parFlankLength+1;
                my $contig_e=$mapinfo{'qs'}+$parFlankLength-1;
                my ($isPassed,$BP_chra,$BP_sitea,$BP_chrb,$BP_siteb)=&Verify_SV($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,$contigID,$contig_s,$contig_e,"CTX");
                if($isPassed == 1){
                    #print "$BP_chra\t$BP_sitea\t$BP_chrb\t$BP_siteb\tCTX\t$contigID\t$last_map{'qe'}\t$mapinfo{'qs'}\n";
                    print $fh_logfile "Approved\n";
                    &Format_Breakpoints($BP_chra,$BP_sitea,$BP_chrb,$BP_siteb,"CTX");
                }
            }
            else{
                &Format_Breakpoints($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,"CTX");
                #print "$BP_ChrA\t$BP_SiteA\t$BP_ChrB\t$BP_SiteB\tCTX\t$contigID\t$last_map{'qe'}\t$mapinfo{'qs'}\n";
            }

        }
        else{
            if($last_map{'dr'} ne $mapinfo{'dr'}){ ##inversion
                pop(@maphits);
                if(@maphits >= 2){
                    &Check_Order(\@maphits);
                }
                @maphits=();

                my $BP_ChrA;
                my $BP_SiteA;
                my $BP_ChrB;
                my $BP_SiteB;
                $BP_ChrA=$last_map{'rn'};
                $BP_ChrB=$mapinfo{'rn'};
                if($last_map{'dr'} eq "+"){
                    $BP_SiteA = $last_map{'se'};
                }
                else{
                    $BP_SiteA = $last_map{'ss'};
                }
                if($mapinfo{'dr'} eq "+"){
                    $BP_SiteB = $mapinfo{'ss'};
                }
                else{
                    $BP_SiteB = $mapinfo{'se'};
                }

                if($parIsVerify == 1){
                    my $contig_s=$last_map{'qe'}-$parFlankLength+1;
                    my $contig_e=$mapinfo{'qs'}+$parFlankLength-1;
                    my ($isPassed,$BP_A,$BP_B)=&Verify_SV($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,$contigID,$contig_s,$contig_e,"INV");
                    if($isPassed == 1){
                        print $fh_logfile "Approved\n";
                        &Format_Breakpoints($BP_ChrA,$BP_A,$BP_ChrB,$BP_B,"INV");
                        #print "$BP_ChrA\t$BP_A\t$BP_ChrB\t$BP_B\tINV\t$contigID\t$last_map{'qe'}\t$mapinfo{'qs'}\n";
                    }
                }
                else{
                    &Format_Breakpoints($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,"INV");
                    #print "$BP_ChrA\t$BP_SiteA\t$BP_ChrB\t$BP_SiteB\tINV\t$contigID\t$last_map{'qe'}\t$mapinfo{'qs'}\n";
                }
            }
        }
        %last_map=%mapinfo;
    }
    if(@maphits >= 2){
        &Check_Order(\@maphits);
    }
}

close $fh_logfile;

sub Check_Order{
    my ($refMaps)=shift;
    my @RefMapinfos=@{$refMaps};

    my %hashSites;
    my %hashLinks;
    ##Create a order according to the contig
    my $direct=$RefMapinfos[0]->{'dr'};
    my $BP_Chr=$RefMapinfos[0]->{'rn'};
    my $ContigName=$RefMapinfos[0]->{'qn'};
    my $contig_Start=$RefMapinfos[0]->{'qs'};
    my $contig_End=$RefMapinfos[scalar(@RefMapinfos)-1]->{'qe'};

    if($direct eq "+"){
        for(my $i=0; $i<@RefMapinfos; $i++){
            my %mapinfo=%{$RefMapinfos[$i]};
            $hashSites{$mapinfo{'ss'}}=$i+1;
            $hashLinks{$mapinfo{'ss'}}=$RefMapinfos[$i];
        }
    }
    else{
        my $count=1;
        for(my $i= scalar(@RefMapinfos)-1; $i>=0;$ i--){
            my %mapinfo=%{$RefMapinfos[$i]};
            $hashSites{$mapinfo{'ss'}}=$count;
            $hashLinks{$mapinfo{'ss'}}=$RefMapinfos[$i];
            $count++;
        }
    }

    my @ITXpoints;
    my @OrderedSites=sort {$a<=>$b} keys %hashSites; ##Re-sort the order according to the reference

    for(my $i=0;$i<@OrderedSites-1;$i++){
        my $currorder=$hashSites{$OrderedSites[$i]};
        my $nextorder=$hashSites{$OrderedSites[$i+1]};
        next if($currorder>$nextorder);
        my $orderdist=$nextorder-$currorder;
        if($orderdist > 1){ ##breakpoint;
            my $breakorder_a=$currorder+1;
            my $breakorder_b=$nextorder-1;
            my $break_ss_a=&Search_MapOrder($breakorder_a,\%hashSites);
            my $break_ss_b=&Search_MapOrder($breakorder_b,\%hashSites);

            my %breakpoint_1_a=%{$hashLinks{$OrderedSites[$i]}};
            my %breakpoint_1_b=%{$hashLinks{$break_ss_a}};

            my %breakpoint_2_a=%{$hashLinks{$break_ss_b}};
            my %breakpoint_2_b=%{$hashLinks{$OrderedSites[$i+1]}};

            my $Chr=$breakpoint_1_a{'rn'};

            my $contig_bp_s;
            my $contig_bp_e;
            if($direct eq "+"){
                #print "$Chr\t$breakpoint_1_a{'se'}\t$Chr\t$breakpoint_1_b{'ss'}\tITX\t$breakpoint_1_a{'qn'}\t$breakpoint_1_a{'qe'}\t$breakpoint_1_b{'qs'}\n";
                $contig_bp_s=$breakpoint_1_a{'qe'}-$parFlankLength+1;
                $contig_bp_e=$breakpoint_1_b{'qs'}+$parFlankLength-1;
            }
            else{
                #print "$Chr\t$breakpoint_1_a{'se'}\t$Chr\t$breakpoint_1_b{'ss'}\tITX\t$breakpoint_1_a{'qn'}\t$breakpoint_1_b{'qe'}\t$breakpoint_1_a{'qs'}\n";
                $contig_bp_s=$breakpoint_1_b{'qe'}-$parFlankLength+1;
                $contig_bp_e=$breakpoint_1_a{'qs'}+$parFlankLength-1;
            }

            if($parIsVerify == 1){
                $contig_bp_s = $contig_Start if($contig_bp_s < $contig_Start);
                $contig_bp_e = $contig_End if($contig_bp_e > $contig_End);
                my($isPassed,$BP_A,$BP_B)=&Verify_SV($BP_Chr,$breakpoint_1_a{'se'},$BP_Chr,$breakpoint_1_b{'ss'},$ContigName,$contig_bp_s,$contig_bp_e,"ITX",$direct);
                if($isPassed == 1){
                    print $fh_logfile "Approved\n";
                    &Format_Breakpoints($BP_Chr,$BP_A,$BP_Chr,$BP_B,"ITX");
                    #print "$BP_Chr\t$BP_A\t$BP_Chr\t$BP_B\tITX\t$ContigName\t$contig_bp_s,$contig_bp_e\n";
                }
            }
            else{
                #if($direct eq "+"){
                #    print "$Chr\t$breakpoint_1_a{'se'}\t$Chr\t$breakpoint_1_b{'ss'}\tITX\t$breakpoint_1_a{'qn'}\t$breakpoint_1_a{'qe'}\t$breakpoint_1_b{'qs'}\n";
                #}
                #else{
                #    print "$Chr\t$breakpoint_1_a{'se'}\t$Chr\t$breakpoint_1_b{'ss'}\tITX\t$breakpoint_1_a{'qn'}\t$breakpoint_1_b{'qe'}\t$breakpoint_1_a{'qs'}\n";
                #}
                &Format_Breakpoints($Chr,$breakpoint_1_a{'se'},$Chr,$breakpoint_1_b{'ss'},"ITX");
            }

            $contig_bp_s="";
            $contig_bp_e="";
            if($direct eq "+"){
                #print "$Chr\t$breakpoint_2_a{'se'}\t$Chr\t$breakpoint_2_b{'ss'}\tITX\t$breakpoint_2_a{'qn'}\t$breakpoint_2_a{'qe'}\t$breakpoint_2_b{'qs'}\n";
                $contig_bp_s=$breakpoint_2_a{'qe'}-$parFlankLength+1;
                $contig_bp_e=$breakpoint_2_b{'qs'}+$parFlankLength-1;
            }
            else{
                #print "$Chr\t$breakpoint_2_a{'se'}\t$Chr\t$breakpoint_2_b{'ss'}\tITX\t$breakpoint_2_a{'qn'}\t$breakpoint_2_b{'qe'}\t$breakpoint_2_a{'qs'}\n";
                $contig_bp_s=$breakpoint_2_b{'qe'}-$parFlankLength+1;
                $contig_bp_e=$breakpoint_2_a{'qs'}+$parFlankLength-1;
            }

            if($parIsVerify == 1){
                $contig_bp_s = $contig_Start if($contig_bp_s < $contig_Start);
                $contig_bp_e = $contig_End if($contig_bp_e > $contig_End);
                my($isPassed,$BP_A,$BP_B)=&Verify_SV($BP_Chr,$breakpoint_2_a{'se'},$BP_Chr,$breakpoint_2_b{'ss'},$ContigName,$contig_bp_s,$contig_bp_e,"ITX",$direct);
                if($isPassed == 1){
                    print $fh_logfile "Approved\n";
                    &Format_Breakpoints($BP_Chr,$BP_A,$BP_Chr,$BP_B,"ITX");
                    #print "$BP_Chr\t$BP_A\t$BP_Chr\t$BP_B\tITX\t$ContigName\t$contig_bp_s,$contig_bp_e\n";
                }
            }
            else{
                #if($direct eq "+"){
                #    print "$Chr\t$breakpoint_2_a{'se'}\t$Chr\t$breakpoint_2_b{'ss'}\tITX\t$breakpoint_2_a{'qn'}\t$breakpoint_2_a{'qe'}\t$breakpoint_2_b{'qs'}\n";
                #}
                #else{
                #    print "$Chr\t$breakpoint_2_a{'se'}\t$Chr\t$breakpoint_2_b{'ss'}\tITX\t$breakpoint_2_a{'qn'}\t$breakpoint_2_b{'qe'}\t$breakpoint_2_a{'qs'}\n";
                #}
                &Format_Breakpoints($Chr,$breakpoint_2_a{'se'},$Chr,$breakpoint_2_b{'ss'},"ITX");
            }
        }
        elsif($orderdist == 1){
            my %curr_map=%{$hashLinks{$OrderedSites[$i]}};
            my %next_map=%{$hashLinks{$OrderedSites[$i+1]}};
            my $refgaplength=$next_map{'ss'}-$curr_map{'se'}+1;
            my $querygaplength;
            if($direct eq "+"){
                $querygaplength=$next_map{'qs'}-$curr_map{'qe'}+1;
            }
            else{
                $querygaplength=$curr_map{'qs'}-$next_map{'qe'}+1;
            }
            if($refgaplength > $parITXLength && $refgaplength/$querygaplength > 5){

                my $contig_bp_s;
                my $contig_bp_e;

                if($direct eq "+"){
                    #print "$curr_map{'rn'}\t$curr_map{'se'}\t$curr_map{'rn'}\t$next_map{'ss'}\tITX\t$curr_map{'qn'}\t$curr_map{'qe'}\t$next_map{'qs'}\n";
                    $contig_bp_s=$curr_map{'qe'}-$parFlankLength+1;
                    $contig_bp_e=$next_map{'qs'}+$parFlankLength-1;
                }
                else{
                    #print "$curr_map{'rn'}\t$curr_map{'se'}\t$curr_map{'rn'}\t$next_map{'ss'}\tITX\t$curr_map{'qn'}\t$next_map{'qe'}\t$curr_map{'qs'}\n";
                    $contig_bp_s=$next_map{'qe'}-$parFlankLength+1;
                    $contig_bp_e=$curr_map{'qs'}+$parFlankLength-1;
                }

                if($parIsVerify == 1){
                    $contig_bp_s = $contig_Start if($contig_bp_s < $contig_Start);
                    $contig_bp_e = $contig_End if($contig_bp_e > $contig_End);
                    my($isPassed,$BP_A,$BP_B)=&Verify_SV($BP_Chr,$curr_map{'se'},$BP_Chr,$next_map{'ss'},$ContigName,$contig_bp_s,$contig_bp_e,"ITX",$direct);
                    if($isPassed == 1){
                        print $fh_logfile "Approved\n";
                        &Format_Breakpoints($BP_Chr,$BP_A,$BP_Chr,$BP_B,"ITX");
                        #print "$BP_Chr\t$BP_A\t$BP_Chr\t$BP_B\tITX\t$ContigName\t$contig_bp_s,$contig_bp_e\n";
                    }
                }
                else{
                    #if($direct eq "+"){
                    #    print "$curr_map{'rn'}\t$curr_map{'se'}\t$curr_map{'rn'}\t$next_map{'ss'}\tITX\t$curr_map{'qn'}\t$curr_map{'qe'}\t$next_map{'qs'}\n";
                    #}
                    #else{
                    #    print "$curr_map{'rn'}\t$curr_map{'se'}\t$curr_map{'rn'}\t$next_map{'ss'}\tITX\t$curr_map{'qn'}\t$next_map{'qe'}\t$curr_map{'qs'}\n";
                    #}
                    &Format_Breakpoints($curr_map{'rn'},$curr_map{'se'},$curr_map{'rn'},$next_map{'ss'},"ITX");
                }
            }
        }
    }

}

sub Search_MapOrder{
    my ($ordernum,$refhashSites)=@_;
    my %hashSites=%{$refhashSites};
    foreach my $key(keys %hashSites){
        return $key if($hashSites{$key} == $ordernum);
    }
}

sub CheckCentro{
    my $refMapinfo=shift;
    my %mapinfo=%{$refMapinfo};
    my $map_rn=$mapinfo{'rn'};
    my $map_ss=$mapinfo{'ss'};
    my $map_se=$mapinfo{'se'};
    my $map_ml=$mapinfo{'ml'};
    my $Centro_S=$CentroRegion{$map_rn}->{'s'};
    my $Centro_E=$CentroRegion{$map_rn}->{'e'};
    if($map_ss >= $Centro_E || $map_se <= $Centro_S){
        return 0;
    }
    elsif($map_se > $Centro_S && $map_se <= $Centro_E){
        my $overlength;
        if($map_ss < $Centro_S){
            $overlength = $map_se - $Centro_S + 1;
            if($overlength/$map_ml >= $parOverpercent){
                return 1;
            }
        }
        else{
            return 1;
        }
    }
    elsif($map_ss >= $Centro_S && $map_ss < $Centro_E){
        my $overlength;
        if($map_se >= $Centro_E){
            $overlength= $Centro_E - $map_ss + 1;
            if($overlength/$map_ml >= $parOverpercent){
                return 1;
            }
        }
        else{
            return 1;
        }
    }
    return 0;
}

sub CheckRepeat{
    my $refMapinfo=shift;
    my %mapinfo=%{$refMapinfo};
    my $map_rn=$mapinfo{'rn'};
    my $map_ss=$mapinfo{'ss'};
    my $map_se=$mapinfo{'se'};
    my $map_ml=$mapinfo{'ml'};
    my @repeats=@{$RepeatRegion{$map_rn}};
    for(my $i=0;$i<@repeats;$i++){
        my %repinfo=%{$repeats[$i]};
        if($map_ss >= $repinfo{'rs'} && $map_ss <= $repinfo{'re'}){
            my $overlength;
            if($map_se >= $repinfo{'re'}){
                $overlength=$repinfo{'re'}-$map_ss+1;
                if($overlength/$map_ml >= $parOverpercent){
                    return 1;
                }
            }
            else{
                return 1;
            }
        }
        elsif($map_se >= $repinfo{'rs'} && $map_se <= $repinfo{'re'}){
            my $overlength;
            if($map_ss <= $repinfo{'rs'}){
                $overlength=$map_se - $repinfo{'rs'} + 1;
                if($overlength/$map_ml >= $parOverpercent){
                    return 1;
                }
            }
            else{
                return 1;
            }
        }
        if($map_se < $repinfo{'re'}){
            return 0;
        }
    }
    return 0;
}

sub MergeSite{
    my @sites=@_;
    my @merged_sites;
    my $last_S=0;
    my $last_E=0;
    my $last_C=0;
    foreach my $refsite (@sites){
        my %siteinfo=%{$refsite};
        if($last_S == 0){
            $last_S = $siteinfo{'rs'};
            $last_E = $siteinfo{'re'};
            $last_C = $siteinfo{'rc'};
            next;
        }
        if($siteinfo{'rs'} >= $last_S && $siteinfo{'rs'} <= $last_E){ ##Overlap
            if($siteinfo{'re'} >= $last_E){
                $last_E=$siteinfo{'re'};
            }
            $last_C+=$siteinfo{'rc'};
        }
        else{
            my %RE_Region;
            $RE_Region{'rs'}=$last_S;
            $RE_Region{'re'}=$last_E;
            $RE_Region{'rc'}=$last_C;
            push(@merged_sites,\%RE_Region);
            $last_S = $siteinfo{'rs'};
            $last_E = $siteinfo{'re'};
            $last_C = $siteinfo{'rc'};
        }
    }
    my %RE_Region;
    $RE_Region{'rs'}=$last_S;
    $RE_Region{'re'}=$last_E;
    $RE_Region{'rc'}=$last_C;
    push(@merged_sites,\%RE_Region);
    return (@merged_sites);
}

sub SearchContig{
    my ($refmapinfo)=shift;
    my %mapinfo=%{$refmapinfo};
    my $map_rn=$mapinfo{'rn'};
    my $map_qs=$mapinfo{'qs'};
    my $map_qe=$mapinfo{'qe'};
    my $contig_Name=$mapinfo{'qn'};
    my @RefMapinfos=@{$RefCoords{$map_rn}};
    foreach my $info (@RefMapinfos){
        my %refinfo=%{$info};
        my $ref_qn=$refinfo{'qn'};
        my $ref_qs=$refinfo{'qs'};
        my $ref_qe=$refinfo{'qe'};
        return 1 if(($contig_Name eq $ref_qn) && ($ref_qs == $map_qs) && ($ref_qe == $map_qe));
    }
    return 0;
}
sub Verify_SV{
    my ($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,$ContigName,$Contig_Start,$Contig_End,$SV_Type,$map_dir)=@_;
    print $fh_logfile "$BP_ChrA\t$BP_SiteA\t$BP_ChrB\t$BP_SiteB\t$ContigName\t$Contig_Start\t$Contig_End\t$SV_Type\n";

    my $seq=$QuerySeqs{$ContigName};
    $Contig_Start=1 if($Contig_Start <= 0);
    my $seq_length=$Contig_End-$Contig_Start;
    my $subseq=substr($seq,$Contig_Start-1,$seq_length);

    my $FasTempfile=File::Temp::tempnam(".","fastmp");
    $FasTempfile=basename($FasTempfile);
    $FasTempfile.=".fas";
    open (my $fh_fasfile,">$FasTempfile");
    print $fh_fasfile ">$ContigName\_$Contig_Start\_$Contig_End\n$subseq\n";
    close $fh_fasfile;
    `$CROSSMATCH $FasTempfile $parRefFile > $FasTempfile.aln`;

    open(my $fh_alnfile,"$FasTempfile.aln");
    my $isReadStart=0;
    my @alns;
    ##re-format cross_over output
    while(<$fh_alnfile>){
        chomp();
	next if($_ eq "");
	if(/^Maximal single base matches/){
            $isReadStart=1;
            next;
        }
	last if(/^\d+\smatching entries/);

        if($isReadStart==1){
            my @lines=split(/\s+/,$_);
            shift(@lines) if($lines[0] eq "");
            if($lines[8] eq "C"){
                $lines[8]="-";
            }
            else{
                splice(@lines,8,0,"+");
            }
            my $newline=join("\t",@lines);
            print $fh_logfile "$newline\n";
            push(@alns,$newline);
        }
    }
    close $fh_alnfile;
    unlink ($FasTempfile,"$FasTempfile.aln","$FasTempfile.log");

    ##------
    if($SV_Type eq "CTX"){
        return 0 if(@alns < 2);
        my ($isPassed,$New_BP_ChromA,$New_BP_SiteA,$New_BP_ChromB,$New_BP_SiteB) = &Check_CTX($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,\@alns);
        if($isPassed == 1){
            return ($isPassed,$New_BP_ChromA,$New_BP_SiteA,$New_BP_ChromB,$New_BP_SiteB);
        }
        else{
            return (0);
        }

    }
    elsif($SV_Type eq "INV"){
        return 0 if(@alns < 2);
        my ($isPassed, $BP_siteA, $BP_siteB) = &Check_INV($BP_ChrA,$BP_SiteA,$BP_SiteB,\@alns);
        if($isPassed == 1){
            return ($isPassed,$BP_siteA,$BP_siteB);
        }
        else{
            return (0);
        }
    }
    else{ ##ITX
        return 0 if(@alns < 2);
        my ($isPassed,$BP_siteA,$BP_siteB) = &Check_ITX($BP_ChrA,$BP_SiteA,$BP_SiteB,\@alns,$map_dir);
        if($isPassed == 1){
            return ($isPassed, $BP_siteA, $BP_siteB);
        }
        else{
            return (0);
        }

    }
}

sub Check_CTX{
    my ($BP_ChrA,$BP_SiteA,$BP_ChrB,$BP_SiteB,$ref_alns)=@_;
    my @alns=@{$ref_alns};
    my $lastChrom;
    my $lastEnd;
    my @breakpoints;

    for(my $i=0;$i<@alns;$i++){
        my @lines=split(/\t/,$alns[$i]);
        my $refN=$lines[9];
        next if($refN ne $BP_ChrA && $refN ne $BP_ChrB);
        next if ($lines[6] - $lines[5] < $parCrossOverMinLength);
        next if ($lines[1] > $parCrossOverMaxSub);
        if($lastChrom eq ""){
            $lastChrom = $refN;
            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
            next;
        }
        if($refN eq $lastChrom){ ##same as before, $refN is BP_ChrA or BP_ChrB
            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
        }
        elsif($refN eq $BP_ChrA || $refN eq $BP_ChrB){ ##different from last
            my %points;
            $points{'chr1'}=$lastChrom;
            $points{'site1'}=$lastEnd;
            $points{'chr2'}=$lines[9];
            if($lines[8] eq "+"){
                $points{'site2'}=$lines[10];
            }
            else{
                $points{'site2'}=$lines[11];
            }
            push(@breakpoints,\%points);

            $lastChrom = $lines[9];

            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
        }
    }

    if(@breakpoints < 1){
        return 0;
    }
    elsif(@breakpoints == 1){
        my %points=%{$breakpoints[0]};
        return (1,$points{'chr1'},$points{'site1'},$points{'chr2'},$points{'site2'});
    }
    else{
        my $pointrec="";
        my $lastdist="";
        for(my $i=0;$i<@breakpoints;$i++){
            my %points=%{$breakpoints[$i]};
            my ($multidist,$dist1,$dist2);
            if($points{'chr1'} eq $BP_ChrA){
                $dist1=abs($points{'site1'}-$BP_SiteA);
                $dist2=abs($points{'site2'}-$BP_SiteB);
                $multidist=$dist1*$dist2;
            }
            else{
                $dist1=abs($points{'site2'}-$BP_SiteA);
                $dist2=abs($points{'site1'}-$BP_SiteB);
                $multidist=$dist1*$dist2;
            }
            if($pointrec eq "" || $multidist < $lastdist){
                $pointrec=$i;
                $lastdist=$multidist;
            }
        }
        my %points=%{$breakpoints[$pointrec]};
        return (1,$points{'chr1'},$points{'site1'},$points{'chr2'},$points{'site2'});
    }

}

sub Check_INV{
    my ($BP_Chr,$BP_SiteA,$BP_SiteB,$ref_alns)=@_;
    my @alns=@{$ref_alns};
    my $last_dir;
    my $lastEnd;
    my @breakpoints;

    if($BP_SiteA > $BP_SiteB){ ##make site A and B in order
        my $temp=$BP_SiteA;
        $BP_SiteA=$BP_SiteB;
        $BP_SiteB=$temp;
    }

    for(my $i=0;$i<@alns;$i++){
        my @lines=split(/\t/,$alns[$i]);
        my $refN=$lines[9];
        next if ($refN ne $BP_Chr);
        next if ($lines[6] - $lines[5] < $parCrossOverMinLength);
        next if ($lines[1] > $parCrossOverMaxSub);
        if($last_dir eq ""){
            $last_dir = $lines[8];
            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
            next;
        }
        if($lines[8] ne $last_dir){
            my %points;
            $points{'site1'}=$lastEnd;
            if($lines[8] eq "+"){
                $points{'site2'}=$lines[10];
            }
            else{
                $points{'site2'}=$lines[11];
            }
            push(@breakpoints,\%points);

            $last_dir = $lines[8];

            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
        }
        else{
            if($lines[8] eq "+"){
                $lastEnd=$lines[11];
            }
            else{
                $lastEnd=$lines[12];
            }
        }
    }

    if(@breakpoints < 1){
        return 0;
    }
    elsif(@breakpoints == 1){
        my %points=%{$breakpoints[0]};
        return (1,$points{'site1'},$points{'site2'});
    }
    else{
        my $pointrec="";
        my $lastdist="";
        for(my $i=0;$i<@breakpoints;$i++){
            my %points=%{$breakpoints[$i]};
            my ($multidist,$dist1,$dist2);
            if($points{'site1'} < $points{'$site2'}){
                $dist1=abs($BP_SiteA-$points{'site1'});
                $dist2=abs($BP_SiteB-$points{'site2'});
                $multidist=$dist1*$dist2;
            }
            else{
                $dist1=abs($BP_SiteA-$points{'site2'});
                $dist2=abs($BP_SiteB-$points{'site1'});
                $multidist=$dist1*$dist2;
            }
            if($pointrec eq "" || $multidist < $lastdist){
                $pointrec=$i;
                $lastdist=$multidist;
            }
        }
        my %points=%{$breakpoints[$pointrec]};
        return (1,$points{'site1'},$points{'site2'});
    }
}

sub Check_ITX{
    my ($BP_Chr,$BP_SiteA,$BP_SiteB,$ref_alns,$map_dir)=@_;
    my @alns=@{$ref_alns};
    my @breakpoints;
    my $lastEnd="";

    if($BP_SiteA > $BP_SiteB){ ##make site A and B in order
        my $temp=$BP_SiteA;
        $BP_SiteA=$BP_SiteB;
        $BP_SiteB=$temp;
    }

    for(my $i=0;$i<@alns;$i++){
        my @lines=split(/\t/,$alns[$i]);
        next if ($lines[8] ne $map_dir || $lines[9] ne $BP_Chr);
        next if ($lines[6] - $lines[5] < $parCrossOverMinLength);
        next if ($lines[1] > $parCrossOverMaxSub);
        if($lastEnd eq ""){
            $lastEnd = $lines[11];
            next;
        }
        else{
            my %points;
            $points{'site1'}=$lastEnd;
            if($map_dir eq "+"){
                $points{'site2'}=$lines[10];
            }
            else{
                $points{'site2'}=$lines[12];
            }
            push (@breakpoints,\%points);
        }
    }

    my $pointrec="";
    my $pointdist="";
    if(@breakpoints < 1){
        return 0;
    }
    else{
        for(my $i=0;$i<@breakpoints;$i++){
            my %breakinfo=%{$breakpoints[$i]};
            my ($distA,$distB,$disttotal);
            if($map_dir eq "+"){
                $distA = abs($breakinfo{'site1'} - $BP_SiteA);
                $distB = abs($breakinfo{'site2'} - $BP_SiteB);
            }
            else{
                $distA = abs($breakinfo{'site2'} - $BP_SiteA);
                $distB = abs($breakinfo{'site1'} - $BP_SiteB);
            }
            if($distA <= $parFlankLength || $distB <= $parFlankLength){
                if($pointrec eq ""){
                    $pointrec=$i;
                    $pointdist=$distA*$distB;
                }
                else{
                    if($pointdist > ($distA*$distB)){
                        $pointrec=$i;
                        $pointdist=$distA*$distB;
                    }
                }
            }
        }

        if($pointrec eq ""){
            return 0;
        }
        else{
            my %point=%{$breakpoints[$pointrec]};
            return (1,$point{'site1'},$point{'site2'});
        }
    }
}

sub Format_Breakpoints{
    my($ChrA,$SiteA,$ChrB,$SiteB,$type)=@_;
    my ($chra,$sitea,$chrb,$siteb);
    
    my $isswap=0;
 
    if(length($ChrA) > length($ChrB)){
        $isswap=1;
    }
    elsif(length($ChrA) == length($ChrB)){
        if($ChrA gt $ChrB){
            $isswap=1;
        }
        elsif($ChrA eq $ChrB){
            if($SiteA > $SiteB){
                $isswap=1;
            }
        }
    }
    
    if($isswap==1){
        $chra=$ChrB;
        $sitea=$SiteB;
        $chrb=$ChrA;
        $siteb=$SiteA;
    }
    else{
        $chra=$ChrA;
        $sitea=$SiteA;
        $chrb=$ChrB;
        $siteb=$SiteB;
    }
    
    if(exists($hashBreakpoints{$chra})){
        my %points=%{$hashBreakpoints{$chra}};
        foreach my $site(sort {$a<=>$b} keys %points){ ##check all the breakpoints
            if(abs($site-$sitea) <= $parMaxBreakpointsOverlap){
                my @breakinfos=@{$points{$site}};
                foreach my $break(@breakinfos){
                    my %breakinfo=%{$break};
                    if($breakinfo{'type'} eq $type){
                        if($breakinfo{'chr'} eq $chrb){
                            if(abs($breakinfo{'site'}-$siteb) <= $parMaxBreakpointsOverlap){
                                return;
                            }
                        }
                    }
                }
            }
        }
        
        my @breakinfos;
        if(exists($points{$sitea})){
            my @breakinfos=@{$points{$sitea}};            
        }
        my %breakinfo;
        $breakinfo{'type'}=$type;
        $breakinfo{'chr'}=$chrb;
        $breakinfo{'site'}=$siteb;
        push(@breakinfos,\%breakinfo);
        $points{$sitea}=\@breakinfos;
        $hashBreakpoints{$chra}=\%points;
    }
    else{
        my @breakinfos;
        my %breakinfo;
        my %points;
        $breakinfo{'type'}=$type;
        $breakinfo{'chr'}=$chrb;
        $breakinfo{'site'}=$siteb;
        push(@breakinfos,\%breakinfo);
        $points{$sitea}=\@breakinfos;
        $hashBreakpoints{$chra}=\%points;
    }
}

#sub Check_ITX{
#    my ($BP_Chr,$ref_alns,$map_dir)=@_;
#    my @alns=@{$ref_alns};
#    my %hashSites;
#    my %hashLinks;
#
#    if($map_dir eq "+"){
#        for(my $i=0;$i<@alns;$i++){
#            my @lines=split(/\t/,$alns[$i]);
#            next if($lines[9]  ne $BP_Chr);
#            next if($lines[8]  ne $map_dir);
#            $hashSites{$lines[10]}=$i+1;
#            $hashLinks{$lines[10]}=$alns[$i];
#        }
#    }
#    else{
#        my $count=1;
#        for(my $i=scalar(@alns)-1;$i>=0;$i--){
#            my @lines=split(/\t/,$alns[$i]);
#            next if($lines[9]  ne $BP_Chr);
#            next if($lines[8]  ne $map_dir);
#            $hashSites{$lines[12]}=$count;
#            $hashLinks{$lines[12]}=$alns[$i];
#            $count++;
#        }
#    }
#
#    my @ITXpoints;
#    my @OrderedSites=sort {$a<=>$b} keys %hashSites; ##Re-sort the order according to the reference
#    my $totalhits=scalar(@OrderedSites)-1;
#    my $isBreakpoints=0;
#    for(my $i=0;$i<@OrderedSites-1;$i++){
#        my $currorder=$hashSites{$OrderedSites[$i]};
#        my $nextorder=$hashSites{$OrderedSites[$i+1]};
#        next if($currorder>$nextorder);
#        my $orderdist=$nextorder-$currorder;
#        if($orderdist > 1){ ##breakpoint;
#            my $breakorder_a=$currorder+1;
#            my $breakorder_b=$nextorder-1;
#            my $break_ss_a=&Search_MapOrder($breakorder_a,\%hashSites);
#            my $break_ss_b=&Search_MapOrder($breakorder_b,\%hashSites);
#
#            my @breakpoint_1_a=split(/\t/,$hashLinks{$OrderedSites[$i]});
#            my @breakpoint_1_b=split(/\t/,$hashLinks{$break_ss_a});
#
#            my @breakpoint_2_a=split(/\t/,$hashLinks{$break_ss_b});
#            my @breakpoint_2_b=split(/\t/,$hashLinks{$OrderedSites[$i+1]});
#
#            my ($BP_a,$BP_b);
#            if($map_dir eq "+"){
#                $BP_a=$breakpoint_1_a[11];
#                $BP_b=$breakpoint_1_b[10];
#            }
#            else{
#                $BP_a=$breakpoint_1_a[12];
#                $BP_b=$breakpoint_1_b[11];
#            }
#            if($parCentroFilter == 1){
#                my $isCentro=&CheckCentro($BP_Chr,$BP_a,$BP_Chr,$BP_b);
#                if($isCentro == 0){
#                    my %points;
#                    $points{'chrom'}=$BP_Chr;
#                    $points{'bp_1'}=$BP_a;
#                    $points{'bp_2'}=$BP_b;
#                    $isBreakpoints=1;
#                    push(@ITXpoints,\%points);
#                }
#            }
#            else{
#                my %points;
#                $points{'chrom'}=$BP_Chr;
#                $points{'bp_1'}=$BP_a;
#                $points{'bp_2'}=$BP_b;
#                $isBreakpoints=1;
#                push(@ITXpoints,\%points);
#            }
#
#            if($map_dir eq "+"){
#                $BP_a=$breakpoint_2_a[11];
#                $BP_b=$breakpoint_2_b[10];
#            }
#            else{
#                $BP_a=$breakpoint_2_a[12];
#                $BP_b=$breakpoint_2_b[11];
#            }
#            if($parCentroFilter == 1){
#                my $isCentro=&CheckCentro($BP_Chr,$BP_a,$BP_Chr,$BP_b);
#                if($isCentro == 0){
#                    my %points;
#                    $points{'chrom'}=$BP_Chr;
#                    $points{'bp_1'}=$BP_a;
#                    $points{'bp_2'}=$BP_b;
#                    $isBreakpoints=1;
#                    push(@ITXpoints,\%points);
#                }
#            }
#            else{
#                my %points;
#                $points{'chrom'}=$BP_Chr;
#                $points{'bp_1'}=$BP_a;
#                $points{'bp_2'}=$BP_b;
#                $isBreakpoints=1;
#                push(@ITXpoints,\%points);
#            }
#        }
#        elsif($orderdist == 1){
#            my $curr_map_info=$hashLinks{$OrderedSites[$i]};
#            my $next_map_info=$hashLinks{$OrderedSites[$i+1]};
#            my $refgaplength;
#            my $querygaplength;
#            my @curr_map=split(/\t/,$curr_map_info);
#            my @next_map=split(/\t/,$next_map_info);
#            my ($BP_a,$BP_b);
#            if($map_dir eq "+"){
#                $BP_a=$curr_map[11];
#                $BP_b=$next_map[10];
#                $refgaplength=$next_map[10]-$curr_map[11]+1;
#                $querygaplength=$next_map[5]-$curr_map[6]+1;
#            }
#            else{
#                $BP_a=$curr_map[12];
#                $BP_b=$next_map[11];
#                $refgaplength=$next_map[11]-$curr_map[12]+1;
#                $querygaplength=$curr_map[5]-$next_map[6]+1;
#            }
#            if($refgaplength > $parITXLength && $refgaplength/$querygaplength > 5){
#                if($parCentroFilter == 1){
#                    my $isCentro=&CheckCentro($BP_Chr,$BP_a,$BP_Chr,$BP_b);
#                    if($isCentro == 0){
#                        my %points;
#                        $points{'chrom'}=$BP_Chr;
#                        $points{'bp_1'}=$BP_a;
#                        $points{'bp_2'}=$BP_b;
#                        $isBreakpoints=1;
#                        push(@ITXpoints,\%points);
#                    }
#                }
#                else{
#                    my %points;
#                    $points{'chrom'}=$BP_Chr;
#                    $points{'bp_1'}=$BP_a;
#                    $points{'bp_2'}=$BP_b;
#                    $isBreakpoints=1;
#                    push(@ITXpoints,\%points);
#                }
#            }
#        }
#    }
#    if($isBreakpoints==1){
#        return (1,\@ITXpoints);
#    }
#    else{
#        return 0;
#    }
#}

sub Usage {#help subprogram
    print << "    Usage";

	Usage: -q <Query_File> -r <Reference_File> -o <Output_File> or -d <Delta_File> -o <Output_File> [options]

        Options: -c     Centromere information

                 -i     Minimal identity requirement (default: 95)

                 -l     Minimal mapping length (default: 500)

                 -p     Minimal coverage ratio (default: 0.8)

                 -t     Minimal length for ITX gap (default: 10000)

                 -f     Flanking length from breakpoint (default: 1000)

                 -v     Breakpoints are verified by cross_over (0:off / 1:on default: off)
                        Must have Query_File and Reference_File to use this function

    Usage

    exit(0);
};
