#!/usr/bin/perl
#===============================================================================
#
#         FILE: Mummer_Tiling.pl
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
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.1
#      CREATED: 22/05/12
#     REVISION: 20/06/12
#===============================================================================

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"c=s","p:s","s:i","r:s","q:s","n:i","help");



if((!defined $opts{c})){
        &Usage();
}


my $corfile=$opts{c};
my $createpsudo=$opts{s};
my $reffile=$opts{r};
my $queryfile=$opts{q};
my $outprefix=$opts{p};
my $isoutputN=$opts{n};

$isoutputN||=0;
$outprefix||=$corfile;

my $tilingfile="$outprefix.tiling";
my $psudofasfile="$outprefix.tiling.fas";

$createpsudo||=0;
if($createpsudo==1){
    &Usage if($queryfile eq "");
    &Usage if($reffile eq "" && $isoutputN==0);
}


my $mincoverage=70;
#my $maxgapratio=5;
my $minindelsize=100;
my $minrevsize=100;
my $minconvsize=100;

my $lastrid;
my $lastreflen;
my @mapcontigs;
my %mapinfos;
open(my $fh_corfile,$corfile) || die "$!\n";
open(my $fh_tilingout,">$tilingfile");

while(<$fh_corfile>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lastrid ne $lines[11]){
        if($lastrid ne ""){
            my $contiginfo=&checkmapping($lastrid,$lastreflen,\@mapcontigs);
        $mapinfos{$lastrid}=$contiginfo;
        }
        $lastrid=$lines[11];
        $lastreflen=$lines[7];
        @mapcontigs=();
    }
    my %mapinfo;
    $mapinfo{rs}=$lines[0];
    $mapinfo{re}=$lines[1];
    $mapinfo{rm}=$lines[4];
    $mapinfo{qm}=$lines[5];
    $mapinfo{mi}=$lines[6];
    $mapinfo{rl}=$lines[7];
    $mapinfo{ql}=$lines[8];
    $mapinfo{rc}=$lines[9];
    $mapinfo{qc}=$lines[10];
    $mapinfo{qn}=$lines[12];
    if($lines[3]>$lines[2]){
        $mapinfo{qd}="+";
        $mapinfo{qs}=$lines[2];
        $mapinfo{qe}=$lines[3];
    }
    else{
        $mapinfo{qd}="-";
        $mapinfo{qs}=$lines[3];
        $mapinfo{qe}=$lines[2];
    }
    push(@mapcontigs,\%mapinfo);
}
my $contiginfo=&checkmapping($lastrid,$lastreflen,\@mapcontigs);
$mapinfos{$lastrid}=$contiginfo;
close $fh_tilingout;

if($createpsudo==1){
    my %refseqs;
    my %queryseqs;
    my $seqid;
    my $seq;

    if($isoutputN==0){
    open(my $fh_reffile,$reffile) || die "$!\n";
    while(<$fh_reffile>){
        chomp();
        if(/^>(\S+)/){
        $refseqs{$seqid}=$seq if($seq ne "");
        $seqid=$1;
        $seq="";
        }
        else{
        $seq.=$_;
        }
    }
    close $fh_reffile;
    $refseqs{$seqid}=$seq;
    }
    $seqid="";
    $seq="";

    open(my $fh_queryfile,$queryfile) || die "$!\n";
    while(<$fh_queryfile>){
    chomp();
    if(/^>(\S+)/){
        $queryseqs{$seqid}=$seq if($seq ne "");
        $seqid=$1;
        $seq="";
    }
    else{
        $seq.=$_;
    }
    }
    close $fh_queryfile;
    $queryseqs{$seqid}=$seq;

    open(my $fh_fasout,">$psudofasfile");

    #foreach my $chrom (sort {$mapinfos{$a}<=>$mapinfos{$b}} keys %mapinfos){
    foreach my $chrom (sort keys %mapinfos){
    my @refmapinfo=@{$mapinfos{$chrom}};
    my @newchrom;
    my $refstart=1;
    my $refseq;
    print ">$chrom\n";
    my $reflength=0;

    if($isoutputN==0){
        if(!exists($refseqs{$chrom})){
            unlink $psudofasfile;
            die "Can't find sequence $chrom in reference file\n";
        }
        $refseq=$refseqs{$chrom};
    }
    my $newchromlength=0;

    for(my $j=0;$j<@refmapinfo;$j++){
        my %mapinfo=%{$refmapinfo[$j]};
        my $maprs=$mapinfo{rs};
        my $mapre=$mapinfo{re};
        my $mapqid=$mapinfo{qn};
        $reflength=$mapinfo{rl} if($reflength==0);
        if($refstart<$maprs){
        my $reflength=$maprs-$refstart;
        $newchromlength+=$reflength;
        my $reffragment;
        if($isoutputN==0){
            $reffragment=substr($refseq,$refstart-1,$reflength);
        }
        else{
            $reffragment=&createNs($reflength);
        }
        print "ref: $refstart,$maprs,$newchromlength\n";
        push(@newchrom,$reffragment);
        }
        my $queryseq;
        if(!exists($queryseqs{$mapqid})){
        unlink $psudofasfile;
        die "Can't find sequence $mapqid in query file\n";
        }

        $queryseq=$queryseqs{$mapqid};

        ##--deal with overlap mapping---
        my $mapqs=1;
        my $mapqe=$mapinfo{ql};
        if($j>0){
        my %lastmapinfo=%{$refmapinfo[$j-1]};
        if($lastmapinfo{re} >= $mapinfo{rs}){
            if($lastmapinfo{mi} >= $mapinfo{mi}){
            my $overlaplength=$lastmapinfo{re}-$mapinfo{rs}+1;
            $mapqs+=$overlaplength;
            }
        }
        }
        if($j<@refmapinfo-1){
        my %nextmapinfo=%{$refmapinfo[$j+1]};
        if($nextmapinfo{rs} <= $mapinfo{re}){
            if($nextmapinfo{mi} > $mapinfo{mi}){
            my $overlaplength=$mapinfo{re}-$nextmapinfo{rs}+1;
            $mapqe-=$overlaplength;
            $refstart=$nextmapinfo{re}+1;
            }
        }
        else{
            $refstart=$mapinfo{re}+1;
        }
        }
        $refstart=$mapinfo{re}+1 if($j==@refmapinfo-1);
        my $querylength=$mapqe-$mapqs+1;
        my $queryfragment=substr($queryseq,$mapqs-1,$querylength);
        $newchromlength+=$querylength;
        print "query: $mapqs,$mapqe,$newchromlength\n";
        if($mapinfo{qd} eq "+"){
        push(@newchrom,$queryfragment);
        }
        else{
        $queryfragment=~tr/atcgATCG/tagcTAGC/;
        $queryfragment= reverse($queryfragment);
        push(@newchrom,$queryfragment);
        }
    }
    my $reffragment;
    if($isoutputN==0){##last reference fragment;
        $reffragment=substr($refseq,$refstart-1);
    }
    else{
        my $length=$reflength-$refstart+1;
        $reffragment=&createNs($length);
    }
    $newchromlength+=length($reffragment);
    print "ref: $refstart,End\n";
    print "Total length: $newchromlength\n";
    push(@newchrom,$reffragment);
    my $newchromseq=join("",@newchrom);
    $newchromseq=&formatseq($newchromseq);
    print $fh_fasout ">$chrom\n$newchromseq";
    }
    close $fh_fasout;
}
exit(1);

sub checkmapping{
    ##steps: 1) check coverage, remove low coverage mapping scaffold.
    ##       2) check linage relationship, keep the same direction info.
    ##       3) check overlap info, solve overlap by length and indentity

    my ($refid,$reflen,$arrayref_mapcontigs)=@_;
    my @contiginfos=@{$arrayref_mapcontigs};
    ##----caculate continous coverage of each contig, remove the contigs' coverage less than mincoverage
    ##----bug may be caused by multi hit in same region
    my $lastcontig;
    my $coverage;
    my @removelines;
    my @currenlines;
    for(my $i=0;$i<@contiginfos;$i++){
        my %mapinfo=%{$contiginfos[$i]};
        if($lastcontig ne $mapinfo{qn}){
            if($coverage>0){
                push (@removelines,@currenlines) if($coverage<$mincoverage);
            }
            $coverage=0;
            @currenlines=();
            $lastcontig=$mapinfo{qn};
        }
        $coverage+=$mapinfo{qc};
        push(@currenlines,$i);
    }
    push (@removelines,@currenlines) if($coverage<$mincoverage);

    my %removelines_hash;
    for(my $i=0;$i<@removelines;$i++){
        $removelines_hash{$removelines[$i]}=1;
    }

    my @fltcontiginfos;
    for(my $i=0;$i<@contiginfos;$i++){
        push(@fltcontiginfos,$contiginfos[$i]) if (!exists($removelines_hash{$i}));
    }
    @contiginfos=@fltcontiginfos;

    my @contigmaps;
    $lastcontig="";
    my %existcontigs;
    my @contigs;

    for(my $i=0;$i<@contiginfos;$i++){
        my %mapinfo=%{$contiginfos[$i]};
        if($lastcontig ne $mapinfo{qn}){
            print "Error: $refid  $mapinfo{qn}\n" if(exists($existcontigs{$mapinfo{qn}}));
        if($lastcontig ne ""){
        my $continfo=&MappingAnalysis(@contigmaps);
        push(@contigs,$continfo);
        }
            @contigmaps=();
            $lastcontig=$mapinfo{qn};
            $existcontigs{$mapinfo{qn}}=1;
        }
        push(@contigmaps,\%mapinfo);
    }
    my $continfo=&MappingAnalysis(@contigmaps);
    push(@contigs,$continfo);

    ##--caculate coverage of reference and new position of the psudo DNA molecular----
    my $totalmaplength;
    my $posshift=0;
    my @newcontigs;
    for(my $i=0;$i<@contigs;$i++){
    my %mapinfo=%{$contigs[$i]};
    $totalmaplength+=$mapinfo{re}-$mapinfo{rs}+1;
    $mapinfo{ns}=$mapinfo{rs}+$posshift;
    my $lengthdiv=$mapinfo{ql}-($mapinfo{re}-$mapinfo{rs}+1);
    $posshift+=$lengthdiv;
    $mapinfo{ne}=$mapinfo{re}+$posshift;
    push(@newcontigs,\%mapinfo);
    }
    my $coverpercent=sprintf("%.2f",$totalmaplength/$reflen*100);
    my $newlength=$reflen+$posshift;
   ##---

   ##---Tiling output---------
    print $fh_tilingout ">$refid $reflen $coverpercent\% $newlength\n";
    for(my $i=0;$i<@newcontigs;$i++){
       my %mapinfo=%{$newcontigs[$i]};
       print $fh_tilingout "$mapinfo{rs}\t$mapinfo{re}\t$mapinfo{ns}\t$mapinfo{ne}\t";
       print $fh_tilingout "$mapinfo{qd}\t$mapinfo{ql}\t$mapinfo{mi}\t$mapinfo{gl}";
       print $fh_tilingout "\t$mapinfo{il}\t$mapinfo{ir}\t$mapinfo{ic}\t$mapinfo{ii}\t$mapinfo{qn}\n";
    }

    return \@newcontigs;
}

sub MappingAnalysis{
    my @contiginfos=@_;
    my %returninfo;
    if(scalar(@contiginfos)==1){
        my %mapinfo=%{$contiginfos[0]};
        $returninfo{rs}=$mapinfo{rs};
        $returninfo{re}=$mapinfo{re};
    $returninfo{rl}=$mapinfo{rl};
    $returninfo{il}=0;
    $returninfo{gl}=0;
        $returninfo{ql}=$mapinfo{ql};
        $returninfo{qn}=$mapinfo{qn};
        $returninfo{qd}=$mapinfo{qd};
        $returninfo{mi}=$mapinfo{mi};
        $returninfo{rm}=$mapinfo{rm};
    $returninfo{ir}=0;
    $returninfo{ic}=0;
    $returninfo{ii}=0;
    }
    else{

        ##calculate main direction of whole mapping
        my $forlength=0;
        my $revlength=0;
        my $qd;
        my $isrev=0;
        for(my $i=0;$i<@contiginfos;$i++){
            my %mapinfo=%{$contiginfos[$i]};
        $qd=$mapinfo{qd} if($qd eq "");
            if($mapinfo{qd} eq "+"){
                $forlength+=$mapinfo{qm};
            }
            else{
                $revlength+=$mapinfo{qm};
            }
        }
        if($forlength > 0 && $revlength > 0){
            if($forlength >= $revlength){
                $isrev=1 if ($revlength > $minrevsize);
                $qd="+";
            }
            else{
                $isrev=1 if ($forlength > $minrevsize);
                $qd="-";
            }
        }

        ##finding fragment reverse

        my $lastctr=0;
        my $isconver=0;
        my $dir1len=0;
        my $dir2len=0;
        for(my $i=0;$i<@contiginfos;$i++){
            my %mapinfo=%{$contiginfos[$i]};
            my $ctrpos=($mapinfo{qs}+$mapinfo{qe})/2;
        if($lastctr==0){
        $lastctr=$ctrpos;
        next;
        }
            if($qd eq "+"){
                if($ctrpos>$lastctr){
                    $lastctr=$ctrpos;
                    if($isconver==0){
                        $dir1len+=$mapinfo{qm};
                    }
                    else{
                        $dir2len+=$mapinfo{qm};
                    }
                }
                else{
                    $isconver=1;
                    $lastctr=$ctrpos;
                }
            }
            else{
                if($ctrpos<$lastctr){
                    $lastctr=$ctrpos;
                    if($isconver==0){
                        $dir1len+=$mapinfo{qm};
                    }
                    else{
                        $dir2len+=$mapinfo{qm};
                    }
                }
                else{
                    $isconver=1;
                    $lastctr=$ctrpos;
                }
            }
        }

        if($dir2len>=$minconvsize && $dir2len>=$minconvsize){
            $isconver=1;
        }
        else{
            $isconver=0;
        }

    my $indelstatus=0; # 1 deletion, 2 insertioin, 3 both
    my $dellength=0;
    my $inslength=0;
        if($isconver==0 && $isrev==0){
        for(my $i=1;$i<@contiginfos;$i++){
        my %mapinfo1=%{$contiginfos[$i-1]};
        my %mapinfo2=%{$contiginfos[$i]};
        my $gaplength=0;
        if($qd eq "+"){
            my $rfspace=$mapinfo2{rs}-$mapinfo1{re};
            $rfspace=0 if($rfspace<0);
            my $qsspace=$mapinfo2{qs}-$mapinfo1{qe};
            $qsspace=0 if($qsspace<0);
            $gaplength=$rfspace-$qsspace;
        }
        else{
            my $rfspace=$mapinfo2{rs}-$mapinfo1{re};
            $rfspace=0 if($rfspace<0);
            my $qsspace=$mapinfo1{rs}-$mapinfo2{re};
            $qsspace=0 if($qsspace<0);
            $gaplength=$rfspace-$qsspace;
        }
        if($gaplength>0){
            $dellength+=$gaplength;
        }
        else{
            $inslength+=abs($gaplength);
        }
        if($gaplength>=$minindelsize){
            if($indelstatus==2){
            $indelstatus=3;
            }
            else{
            $indelstatus=1;
            }
        }
        elsif($gaplength<=(0-$minindelsize)){
            if($indelstatus==1){
            $indelstatus=3;
            }
            else{
            $indelstatus=2;
            }
        }
        }
        }

    my $matchlength;
    my $totallength;
    my $refmatchlength;
    my $lastrfend;

    for(my $i=0;$i<@contiginfos;$i++){
        my %mapinfo=%{$contiginfos[$i]};
        $matchlength+=$mapinfo{mi}*$mapinfo{qm}/100;
        $totallength+=$mapinfo{qm};
        if($lastrfend > $mapinfo{rs}){
        $refmatchlength+=$mapinfo{re}-$lastrfend;
        }
        else{
        $refmatchlength+=$mapinfo{re}-$mapinfo{rs}+1;
        }
        if($returninfo{rs}==0){
        $returninfo{rs}=$mapinfo{rs};
        $returninfo{re}=$mapinfo{re};
        }
        $returninfo{rs}=$mapinfo{rs} if($returninfo{rs} > $mapinfo{rs});
        $returninfo{re}=$mapinfo{re} if($returninfo{re} < $mapinfo{re});
        $lastrfend=$mapinfo{re};
    }
    $returninfo{gl}=$dellength;
    $returninfo{il}=$inslength;
    my %tempinfo=%{$contiginfos[0]};
    $returninfo{ql}=$tempinfo{ql};
    $returninfo{qn}=$tempinfo{qn};
    $returninfo{rl}=$tempinfo{rl};
    $returninfo{qd}=$qd;
    $returninfo{mi}=sprintf("%.2f",$matchlength/$totallength*100);
    $returninfo{rm}=$matchlength;
    $returninfo{ir}=$isrev;
    $returninfo{ic}=$isconver;
    $returninfo{ii}=$indelstatus;
    }
    return \%returninfo;
}

sub createNs{
    my $Nnumber=shift;
    my $Ns;
    for(my $i=0;$i<$Nnumber;$i++){
    $Ns.="N";
    }
    return $Ns;
}

sub formatseq{
    my $seq=shift;
    my $seqlength = length($seq);
    my $modulus = $seqlength % 80;
    $seq =~ s/(.{80})/$1\n/ig;
    $seq="$seq\n" if($modulus!=0);
    return $seq;
}


sub Usage{
    print << "    Usage";

    Usage: $0 <options>

        -c     Results of show-coords (MUMmer)

        -p     Prefixed output file (default: coords file)

        -s     Create psudo molecular by mapping information (default: 0)

        -n     Create psudo molecular fill gaps with N (default:0)

        -r     Reference file (must be given when -n 1 -s 1)

        -q     Query file (must be given when -s 1)

        -help  Show help

    Usage

    exit(0);
}
