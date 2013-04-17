#!/usr/bin/perl
#===============================================================================
#
#         FILE:SV_BreakCheck.pl
#
#        USAGE:SV_BreakCheck.pl -i <Breakpoint File> -o <Output File> --ref <Reference Genome> --bam <Bam File>
#
#  DESCRIPTION:The script is used to verify the breakpoints using denovo assembly and smith-waterman alignment
#
#      OPTIONS: ---
# REQUIREMENTS: Bioperl, velvet, samtools, bam2fastq, cross_over
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 02-20-2013
#     REVISION:
#===============================================================================
use strict;
use Getopt::Long;
use File::Path;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Bio::SearchIO;
use Bio::SeqIO;

my ($BrkFile,$OutFile,$RefFile,$BamFile,$Flankin,$Flankout,$TempDir,$velvet_insert_len,$MaxCoverage);

GetOptions('i=s'=>\$BrkFile,'o=s'=>\$OutFile,'ref=s'=>\$RefFile,
           'bam=s'=>\$BamFile,'cov:i'=>\$MaxCoverage,'in:i'=>\$Flankin,'tmp:s'=>\$TempDir,'out:i'=>\$Flankout,'ins:i'=>\$velvet_insert_len);
if($BrkFile eq "" || $OutFile eq "" || $velvet_insert_len eq "" ||$RefFile eq "" || $BamFile eq ""){
    &Usage();
}

$Flankin ||= 100;
$Flankout ||= 500;
$MaxCoverage ||= 1000;
$velvet_insert_len ||= 200;
$TempDir ||= ".";

my $LogFile="$OutFile.fas";
my $velvet_kmer=31;
my $velvet_cov_cutoff=5;
my $velvet_min_contig=200;
my $Ref_Flanking=$Flankout+200;
my $Contig_min_coverage=60;
my $SV_length_per=10;
my $SV_length_var=$SV_length_per/100;
my $SV_MinBp_length=50;
my $SV_Bp_window=100;
my $Mapping_overlap=0.1;
my $Merge_Flanking=50;

my %Refregion;
my %RefLengths;
my @SVPoints;
my @SVnotpassed;

my $samtools = `which samtools`;
my $bam2fastq = `which bam2fastq`;
my $velveth = `which velveth`;
my $velvetg = `which velvetg`;
my $crossmatch = `which cross_match`;

if($samtools eq ""){
    print "Can't find samtools\n";
    exit(-1);
}
if($bam2fastq eq ""){
    print "Can't find bam2fastq\n";
    exit(-1);
}
if($velveth eq ""){
    print "Can't find velveth\n";
    exit(-1);
}
if($velvetg eq ""){
    print "Can't find velvetg\n";
    exit(-1);
}
if($crossmatch eq ""){
    print "Can't find cross_match\n";
    exit(-1);
}


if(!-e($TempDir)){
    mkdir $TempDir;
}

die "Can't find input file: $BrkFile\n" if (!-e($BrkFile));
die "Can't find bam file: $BamFile\n" if(!-e($BamFile));
die "Can't find reference file: $RefFile\n" if(!-e($RefFile));

my $Bam = Bio::DB::Bam->open($BamFile);
if(!-e("$BamFile.bai")){
    my $status = Bio::DB::Bam->index_build($BamFile);
    if($status < 0){
        print "Can't build index for the bam file: $BamFile\n";
        exit(-1);
    }
}
my $Bam_index = Bio::DB::Bam->index_open($BamFile);
my $Bam_header = $Bam -> header();

{
    my $refnum=$Bam_header -> n_targets();
    for(my $i=0;$i<$refnum;$i++){
        my $refN=$Bam_header->target_name->[$i];
        my $refLength=$Bam_header->target_len->[$i];
        $RefLengths{$refN}=$refLength;
    }
}

open(my $fh_log,">$LogFile");
open(my $fh_brks,"$BrkFile");
while(<$fh_brks>){
    chomp();
    my @lines=split(/\t/,$_);
    my ($ispass,$breakpoints)=checkSV(@lines[0..4]);
    if($ispass==1){
        my ($isoverlap,$index)=&checkBPOverlap($breakpoints);
        if($isoverlap==0){
            $breakpoints->{'soft'}=$lines[5];
            push (@SVPoints,$breakpoints);
        }
        elsif($isoverlap==1){
            $SVPoints[$index]->{'soft'}=$SVPoints[$index]->{'soft'}.":".$lines[5];
        }
    }
    else{
        my %breakinfo;
        $breakinfo{'chrA'}=$lines[0];
        $breakinfo{'chrB'}=$lines[2];
        $breakinfo{'posA'}=$lines[1];
        $breakinfo{'posB'}=$lines[3];
        $breakinfo{'sv'}=$lines[4];
        $breakinfo{'soft'}=$lines[5];
        push(@SVnotpassed,\%breakinfo);
    }
}

close $fh_log;
close $fh_brks;
open(my $fh_out,">$OutFile");

for(my $i=0;$i<@SVPoints;$i++){
    print $fh_out $SVPoints[$i]->{'chrA'},"\t",$SVPoints[$i]->{'posA'},"\t",$SVPoints[$i]->{'chrB'},
    "\t",$SVPoints[$i]->{'posB'},"\t",$SVPoints[$i]->{'sv'},"\tpass\t",$SVPoints[$i]->{'soft'},"\n";
}

for(my $i=0;$i<@SVnotpassed;$i++){
    print $fh_out $SVnotpassed[$i]->{'chrA'},"\t",$SVnotpassed[$i]->{'posA'},"\t",$SVnotpassed[$i]->{'chrB'},
    "\t",$SVnotpassed[$i]->{'posB'},"\t",$SVnotpassed[$i]->{'sv'},"\t.\t",$SVnotpassed[$i]->{'soft'},"\n";
}

close $fh_out;
exit(0);

sub checkBPOverlap{
    my $breakpoints=shift;
    for(my $i=0;$i<@SVPoints;$i++){
        if($SVPoints[$i]->{'sv'} eq $breakpoints->{'sv'}){
            if($SVPoints[$i]->{'chrA'} eq $breakpoints->{'chrA'}
               && $SVPoints[$i]->{'chrB'} eq $breakpoints->{'chrB'}){
                if($breakpoints->{'posA'} >= $SVPoints[$i]->{'posA'}-$Merge_Flanking
                   && $breakpoints->{'posA'} <= $SVPoints[$i]->{'posA'}+$Merge_Flanking){
                    if($breakpoints->{'posA'} >= $SVPoints[$i]->{'posA'}-$Merge_Flanking
                       && $breakpoints->{'posA'} <= $SVPoints[$i]->{'posA'}+$Merge_Flanking){
                        return (1,$i);
                    }
                }
            }
        }
    }
    return 0;
}

sub checkSV{
    my ($Chr1,$Pos1,$Chr2,$Pos2,$SVtype)=@_;

    die "Can't find $Chr1 in bam file: $BamFile\n" if(!exists($RefLengths{$Chr1}));
    die "Can't find $Chr2 in bam file: $BamFile\n" if(!exists($RefLengths{$Chr2}));

    my $Chr1_len=$RefLengths{$Chr1};
    my $Chr2_len=$RefLengths{$Chr2};

    my $Tempfile=File::Temp::tempnam(".","Temp");
    $Tempfile=$TempDir."/$Tempfile";
    my $TempBam=$Tempfile.".bam";

    my $ispass=0;
    my %breakpoints;
    my %velvet_contigs;

    if($SVtype ne "CTX"){
        die "$SVtype: chr2 should be the same as chr1\n" if($Chr1 ne $Chr2);
        if($Pos2 < $Pos1){
            my $postmp=$Pos1;
            $Pos1=$Pos2;
            $Pos2=$postmp;
        }
    }
    my $editedsites;
    my ($Pos1s,$Pos1e,$Pos2s,$Pos2e);
    if($SVtype eq "DEL" ){
        $Pos1s=$Pos1-$Flankout;
        $Pos1s = 1 if($Pos1s < 1);
        $Pos1e=$Pos1+$Flankin;
        $Pos1e = $Chr1_len if($Pos1e > $Chr1_len);
        $Pos2s=$Pos2-$Flankin;
        $Pos2s = 1 if($Pos2s < 1);
        $Pos2e=$Pos2+$Flankout;
        $Pos2e = $Chr2_len if($Pos2e > $Chr2_len);
    }
    elsif($SVtype eq "CTX" || $SVtype eq "INV" ||$SVtype eq "DUP"){
        $Pos1s=$Pos1-$Flankout;
        $Pos1s = 1 if($Pos1s < 1);
        $Pos1e=$Pos1+$Flankout;
        $Pos1e = $Chr1_len if($Pos1e > $Chr1_len);
        $Pos2s=$Pos2-$Flankout;
        $Pos2s = 1 if($Pos2s < 1);
        $Pos2e=$Pos2+$Flankout;
        $Pos2e = $Chr2_len if($Pos2e > $Chr2_len);
    }

    my $isRepReg1=getCoverage($Chr1,$Pos1s,$Pos1e);
    my $isRepReg2=getCoverage($Chr2,$Pos2s,$Pos2e);

    if($isRepReg1 > 0 || $isRepReg2 > 0){
        #print "Repeat Region\n";
        return 0;
    }

    if($SVtype ne "CTX" && $SVtype ne "DUP"){
        if($Pos1e>=$Pos2s){
            $editedsites="$Chr1:$Pos1s\-$Pos2e";
        }
        else{
        $editedsites="$Chr1:$Pos1s\-$Pos1e $Chr2:$Pos2s-$Pos2e";
        }
    }
    elsif($SVtype eq "DUP"){
        $editedsites="$Chr1:$Pos1s\-$Pos1e";
    }
    elsif($SVtype eq "CTX"){
        $editedsites="$Chr1:$Pos1s\-$Pos1e $Chr2:$Pos2s-$Pos2e";
    }

    {
        my $RefFasta;
        my $RefGenome=Bio::DB::Fasta->new($RefFile);
        if($SVtype ne "CTX"){
            my $seqS=$Pos1-$Ref_Flanking;
            $seqS=1 if ($seqS < 1);
            my $seqE=$Pos2+$Ref_Flanking;
            $seqE=$Chr1_len if($seqE>$Chr1_len);
            my %seqinfo={"chr"=>$Chr1,"s"=>$seqS,"e"=>$seqE};

            $Refregion{$Chr1}->{'s'}=$seqS;
            $Refregion{$Chr1}->{'e'}=$seqE;

            my $refregion="$Chr1:$seqS-$seqE";
            my $refseq=$RefGenome->seq($refregion);
            my $refname="$Chr1\_$seqS\-$seqE";
            $RefFasta=">$refname\n$refseq\n";
        }
        else{
            my $seq1s=$Pos1-$Ref_Flanking;
            $seq1s=1 if ($seq1s < 1);
            my $seq1e=$Pos1+$Ref_Flanking;
            $seq1e=$Chr1_len if($seq1e > $Chr1_len);

            my $seq2s=$Pos2-$Ref_Flanking;
            $seq2s=1 if ($seq2s < 1);
            my $seq2e=$Pos2+$Ref_Flanking;
            $seq2e=$Chr2_len if($seq2e > $Chr2_len);

            $Refregion{$Chr1}->{'s'}=$seq1s;
            $Refregion{$Chr1}->{'e'}=$seq1e;
            $Refregion{$Chr2}->{'s'}=$seq2s;
            $Refregion{$Chr2}->{'e'}=$seq2e;

            my $refregion1="$Chr1:$seq1s-$seq1e";
            my $refregion2="$Chr2:$seq2s-$seq2e";

            my $refseq1=$RefGenome->seq($refregion1);
            my $refseq2=$RefGenome->seq($refregion2);

            my $refname1="$Chr1\_$seq1s\-$seq1e";
            my $refname2="$Chr2\_$seq2s\-$seq2e";

            $RefFasta=">$refname1\n$refseq1\n>$refname2\n$refseq2\n";
        }

        open(my $fh_outref,">$Tempfile.ref.fasta");
        print $fh_outref $RefFasta;
        close $fh_outref;
    }

    {
        `samtools view -bh $BamFile $editedsites > $TempBam`;
        `bam2fastq -q -o $Tempfile\#.fastq $TempBam`;
        `velveth $Tempfile $velvet_kmer -fastq -shortPaired $Tempfile\_1.fastq $Tempfile\_2.fastq`;
        `velvetg $Tempfile -cov_cutoff $velvet_cov_cutoff -min_contig_lgth $velvet_min_contig -ins_length $velvet_insert_len -exp_cov auto`;
        my($isFas,$ref_fas)=&readFasta("$Tempfile\/contigs.fa");
        if($isFas == 0){
            goto FINISH;
        }
        %velvet_contigs=%{$ref_fas};
        `cross_match $Tempfile\/contigs.fa $Tempfile.ref.fasta > $Tempfile.aln`;
    }

    my %cross_res=%{&readCrossMatch("$Tempfile.aln",\%velvet_contigs)};
    my $matchcontig;

    if($SVtype eq "DEL" ||$SVtype eq "DUP"){
        for my $contigName (keys %cross_res){
            my @mapinfo=@{$cross_res{$contigName}};
            my $ispassed=0;
            next if (@mapinfo == 1);
            my $mapping_dir=mapDirection(@mapinfo);
            next if ($mapping_dir eq "m");
            my $subName=$mapinfo[0]->{'sN'};
            my $subS=$Refregion{$subName}->{'s'};
            for(my $i=0;$i<@mapinfo-1;$i++){
                my $q0e=$mapinfo[$i]->{'qE'};
                my $q1s=$mapinfo[$i+1]->{'qS'};
                my $overlen=abs($q1s-$q0e);
                my $contiglen=length($velvet_contigs{$contigName});
                next if($overlen/$contiglen > $Mapping_overlap);

                my $pointA=$mapinfo[$i]->{'sE'}+$subS-1; #0E
                my $pointB=$mapinfo[$i+1]->{'sS'}+$subS-1; #1S
                my $breaklength=$Pos2-$Pos1+1;

                if($pointA > $pointB){
                    ($pointA,$pointB)=&switchAB($pointA,$pointB);
                }
                my $bplength=$pointB-$pointA+1;
                if($bplength >= $SV_MinBp_length){
                    if($bplength >= $breaklength*(1-$SV_length_var) && $bplength <= $breaklength*(1+$SV_length_var)){
                        $breakpoints{'sv'}=$SVtype;
                        $breakpoints{'chrA'}=$Chr1;
                        $breakpoints{'chrB'}=$Chr2;
                        $breakpoints{'posA'}=$pointA;
                        $breakpoints{'posB'}=$pointB;
                        $ispass=1;
                        $matchcontig=$contigName;
                        goto FINISH;
                        #print $fh_out "$Chr1\t$pointA\t$Chr2\t$pointB\tpassed\t$Chr1\t$Pos1\t$Chr2\t$Pos2\n";
                    }
                }
            }
        }
    }
    elsif($SVtype eq "INV"){
        ##search for breakpoints, the mapping direction around breakpoints should be different;
        ##so each contigs should have at least two res, in two direction
        my ($brkpos1,$brkpos2);
        for my $contigName (keys %cross_res){
            my @mapinfo=@{$cross_res{$contigName}};
            next if (@mapinfo == 1);
            my $mapping_dir=mapDirection(@mapinfo);
            next if ($mapping_dir ne "m");
            for(my $i=0;$i<@mapinfo-1;$i++){
                my $currdir=$mapinfo[$i]->{'str'};
                my $nextdir=$mapinfo[$i+1]->{'str'};
                my $subName=$mapinfo[0]->{'sN'};
                my $subS=$Refregion{$subName}->{'s'};

                my $q0e=$mapinfo[$i]->{'qE'};
                my $q1s=$mapinfo[$i+1]->{'qS'};
                my $overlen=abs($q1s-$q0e);
                my $contiglen=length($velvet_contigs{$contigName});
                next if($overlen/$contiglen > $Mapping_overlap);

                if($currdir ne $nextdir){  ## 0E - 1S
                    my $pointA=$mapinfo[$i]->{'sE'}+$subS-1; #0E
                    my $pointB=$mapinfo[$i+1]->{'sS'}+$subS-1; #1S
                    my $breaklength=$Pos2-$Pos1+1;
                    if($pointA > $pointB){
                        ($pointA,$pointB)=&switchAB($pointA,$pointB);
                    }
                    my $bplength=$pointB-$pointA+1;
                    my $alignbreaklength=$mapinfo[$i+1]->{'qS'}-$mapinfo[$i]->{'qE'}+1;
                    if($bplength >= $SV_MinBp_length && ($alignbreaklength*5) < $breaklength){
                        if($bplength >= $breaklength*(1-$SV_length_var) && $bplength <= $breaklength*(1+$SV_length_var)){
                            $breakpoints{'sv'}=$SVtype;
                            $breakpoints{'chrA'}=$Chr1;
                            $breakpoints{'chrB'}=$Chr2;
                            $breakpoints{'posA'}=$pointA;
                            $breakpoints{'posB'}=$pointB;
                            $ispass=1;
                            $matchcontig=$contigName;
                            goto FINISH;
                        }
                    }
                }
            }
        }
    }
    elsif($SVtype eq "CTX"){
        for my $contigName (keys %cross_res){
            my @mapinfo=@{$cross_res{$contigName}};
            next if(@mapinfo == 1);
            for(my $i=0;$i<@mapinfo-1;$i++){
                my $subNameA=$mapinfo[$i]->{'sN'};
                my $subSA=$Refregion{$subNameA}->{'s'};
                my $subNameB=$mapinfo[$i+1]->{'sN'};
                my $subSB=$Refregion{$subNameB}->{'s'};
                my $pointA=$mapinfo[$i]->{'sE'}+$subSA-1;
                my $pointB=$mapinfo[$i+1]->{'sS'}+$subSB-1;
                my $q0e=$mapinfo[$i]->{'qE'};
                my $q1s=$mapinfo[$i+1]->{'qS'};
                my $overlen=abs($q1s-$q0e);
                my $contiglen=length($velvet_contigs{$contigName});
                next if($overlen/$contiglen > 0.1);
                if($subNameA ne $subNameB){
                    ($pointA,$pointB)=switchAB($pointA,$pointB) if($subNameA ne $Chr1);
                    $breakpoints{'sv'}=$SVtype;
                    $breakpoints{'chrA'}=$Chr1;
                    $breakpoints{'chrB'}=$Chr2;
                    $breakpoints{'posA'}=$pointA;
                    $breakpoints{'posB'}=$pointB;
                    $ispass=1;
                    $matchcontig=$contigName;
                    goto FINISH;
                }
            }
        }
    }

FINISH:
    rmtree ($Tempfile);
    unlink ("$Tempfile\_1.fastq","$Tempfile\_2.fastq","$Tempfile.ref.fasta","$Tempfile.aln","$Tempfile.bam");
    if($ispass==1){
        print $fh_log ">$Chr1\_$Pos1\_$Chr2\_$Pos2|$Chr1\_$breakpoints{'posA'}\_$Chr2\_$breakpoints{'posB'}|$SVtype\n";
        print $fh_log $velvet_contigs{$matchcontig},"\n";
    }
    return ($ispass,\%breakpoints);
}

sub readCrossMatch{
    my ($alnFile,$ref_fas)=@_;
    my @hits;
    my $lastQueryN;
    my $searchIO = Bio::SearchIO->new( -format => 'cross_match', -file => "$alnFile" );
    my %res;
    while(my $res = $searchIO->next_result) {
        while(my $hit = $res->next_hit) {
            while(my $hsp = $hit->next_hsp) {
                my $queryN=$hsp->{'QUERY_NAME'};
                my @subName=split(/\_/,$hsp->{'HIT_NAME'});
		my $subN=$subName[0];
                my $queryS=$hsp->{'QUERY_START'};
                my $queryE=$hsp->{'QUERY_END'};
                my $convlen=$hsp->{'CONSERVED'};
                my $idpercent=($convlen/($queryE-$queryS))*100;
                my $subS=$hsp->{'HIT_START'};
                my $subE=$hsp->{'HIT_END'};
                my $strand;
                if($subS > $subE){
                    #($subS,$subE)=&switchAB($subS,$subE);
                    $strand='-';
                }
                else{
                    $strand='+';
                }
                if($lastQueryN ne $queryN){
                    if($lastQueryN ne ""){
                        my @sorthits=&hitSort(@hits);
                        my $coverage=&hitCoverpercent($lastQueryN,$ref_fas,@sorthits);
                        if($coverage >= $Contig_min_coverage){
                            $res{$lastQueryN}=\@sorthits;
                        }
                        else{
                            print "$lastQueryN\tCoverage:$coverage\n";
                        }
                        @hits=();
                    }
                    $lastQueryN=$queryN;
                }
                my %hitinfo;
                $hitinfo{'qS'}=$queryS;
                $hitinfo{'qE'}=$queryE;
                $hitinfo{'sN'}=$subN;
                $hitinfo{'sS'}=$subS;
                $hitinfo{'sE'}=$subE;
                $hitinfo{'cL'}=$convlen;
                $hitinfo{'iP'}=$idpercent;
		$hitinfo{'str'}=$strand;
                push(@hits,\%hitinfo);
            }
        }
    }
    my @sorthits=&hitSort(@hits);
    my $coverage=&hitCoverpercent($lastQueryN,$ref_fas,@sorthits);
    if($coverage >= $Contig_min_coverage){
        $res{$lastQueryN}=\@sorthits;
    }
    return \%res;
}

sub mapDirection{
    my @hits=@_;
    my $dir=$hits[0]->{'str'};
    for(my $i=1;$i<@hits;$i++){
	my $hitdir=$hits[$i]->{'str'};
	$dir='m' if($dir ne $hitdir);
    }
    return $dir;
}

sub hitCoverpercent{
    my ($SeqN,$ref_fas,@hits)=@_;
    my %velvet_contigs=%{$ref_fas};
    my $SeqLength=length($velvet_contigs{$SeqN});
    my $lastEnd=0;
    my $coverlength;
    for(my $index=0;$index<@hits;$index++){
        my %hitinfo=%{$hits[$index]};
        if($hitinfo{'qS'} > $lastEnd){
            $coverlength+=$hitinfo{'qE'}-$hitinfo{'qS'}+1;
            $lastEnd = $hitinfo{'qE'};
        }
        else{
            if($hitinfo{'qE'} > $lastEnd){
                $coverlength += $hitinfo{'qE'} - $lastEnd;
                $lastEnd=$hitinfo{'qE'};
            }
        }
    }
    if ($SeqLength == 0) {
        return 0;
    }
    else{
        my $coverpercent= $coverlength/$SeqLength * 100;
        return $coverlength;
    }
}

sub readFasta{
    my $fasFile=shift;
    if(!-e($fasFile)){
        return 0;
    }
    my $seqIn = Bio::SeqIO->new(-file => $fasFile, -format => 'fasta');
    my %fasta;
    my $seqcount=0;
    while(my $seqinfo=$seqIn->next_seq()){
        my $seqID=$seqinfo->display_id;
        my $seq=$seqinfo->seq;
        $fasta{$seqID}=$seq;
        $seqcount++;
    }
    if($seqcount==0){
        return 0;
    }
    else{
        return (1,\%fasta);
    }
}

sub hitSort{
    my @hits=@_;
    my $not_complete = 1;
    my $index;
    return @hits if(@hits<2);
    my $len = ((scalar @hits)-2);
    while($not_complete){
        $not_complete=0;
        foreach $index (0..$len){
            my %hitinfo=%{$hits[$index]};
            my %hitinfo_n=%{$hits[$index+1]};
            if($hitinfo{'qS'} > $hitinfo_n{'qS'}){
                $hits[$index+1]=\%hitinfo;
                $hits[$index]=\%hitinfo_n;
                $not_complete = 1;
            }
        }
    }
    return @hits;
}

sub switchAB{
    my ($numA,$numB)=@_;
    my $numtemp=$numA;
    $numA=$numB;
    $numB=$numtemp;
    return($numA,$numB);
}


sub getCoverage{
    my ($seqid,$start,$end)=@_;
    my $coverage=$Bam_index->coverage($Bam,$seqid,$start,$end);
    my $totalcover;
    for(@{$coverage}){
        $totalcover+=$_;
    }
    my $average_cov=$totalcover/($end-$start+1);
    if($average_cov < $MaxCoverage){
        return 0;
    }
    else{
        return 1;
    }
}


sub Usage #help subprogram
{
    print << "    Usage";

    Usage: $0 -i <Breakpoint File> -o <Output File> --ref <Reference Genome> --bam <Bam File> [options]

        Options: --in <interger>    : Flanking inner length, default 100

                 --out <integer>    : Flanking outer length, default 500

                 --cov <integer>    : Max coverage of the breakpoint flanking region, default 1000

                 --ins <integer>    : Insert length of sequencing library, default 200

                 --tmp <string>     : Temporary directory

    Usage

    exit(0);
}
