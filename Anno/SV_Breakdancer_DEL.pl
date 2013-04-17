#!/usr/bin/perl
use strict;
use File::Basename;
my $maxdel=1000000;
my $mindel=50;
my $minscore=90;
my $minreadssupport=5;
my $mergeoverlap=0.5;

my ($BrkFile,$Sample_Name,$REF_File)=@ARGV;
die "Usage:$0 <Breakdancer_output> <Sample_name> <Reference_name>" if(@ARGV<3);


my $Tempfile1=basename($BrkFile).".tmp1";
my $Tempfile2=basename($BrkFile).".tmp2";
my $mergefile=basename($BrkFile).".DEL.mrg";
my $VCFfile=basename($BrkFile).".DEL.vcf";

open (my $fh_brkfile,"$BrkFile") || die "Can't open file: $BrkFile\n";
open (my $fh_temp1,">$Tempfile1");

while(<$fh_brkfile>){
    next if(/^#/);
    chomp();
    my @lines=split(/\t/,$_);
    next if($lines[6] ne "DEL");
    my $Dellength=$lines[4]-$lines[1];
    next if($Dellength > $maxdel || $Dellength < $mindel || $lines[8] < $minscore);
    print $fh_temp1 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";
}
close $fh_brkfile;
close $fh_temp1;

`cat $Tempfile1 | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $Tempfile2`;


open (my $fh_filein,$Tempfile2) || die "Can't open file $Tempfile2\n";
open (my $fh_mergefile,">$mergefile");
my $lastchrom;
my $lastS;
my $lastE;
my $lastsupportreads;
while(<$fh_filein>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[0] eq $lastchrom){
        if($lines[1] >= $lastE){
            print $fh_mergefile "$lastchrom\t$lastS\t$lastchrom\t$lastE\tDEL\t$lastsupportreads\n";
            $lastchrom=$lines[0];
            $lastS=$lines[1];
            $lastE=$lines[3];
            $lastsupportreads=$lines[7];
        }
        if($lines[3] <= $lastE){
            $lastsupportreads+=$lines[7];
        }
        else{
            my $overlap=$lastE-$lines[1];
            my $lastlength=$lastE-$lastS;
            my $newlength=$lines[3]-$lines[1];
            if($overlap/$lastlength >= $mergeoverlap && $overlap/$newlength >= $mergeoverlap){
                $lastE=$lines[3];
                $lastsupportreads+=$lines[7];
            }
            else{
                print $fh_mergefile "$lastchrom\t$lastS\t$lastchrom\t$lastE\tDEL\t$lastsupportreads\n";
                $lastchrom=$lines[0];
                $lastS=$lines[1];
                $lastE=$lines[3];
                $lastsupportreads=$lines[7];
            }
        }
    }
    else{
        print $fh_mergefile "$lastchrom\t$lastS\t$lastchrom\t$lastE\tDEL\t$lastsupportreads\n" if($lastchrom ne "");
        $lastchrom=$lines[0];
        $lastS=$lines[1];
        $lastE=$lines[3];
        $lastsupportreads=$lines[7];
    }
}
print $fh_mergefile "$lastchrom\t$lastS\t$lastchrom\t$lastE\tDEL\t$lastsupportreads\n";
close $fh_filein;
close $fh_mergefile;

my %RefSeq;

open(my $fh_ref,"$REF_File") ||  die "Can't open file $REF_File";
my $seqN;
my $seqS;
my @SeqNs;
my @SeqLens;
my $reference=basename($REF_File);
while(<$fh_ref>){
    chomp();
    if(/^>(\S+)/){
        if($seqN ne ""){
            $RefSeq{$seqN}=$seqS;
            push(@SeqLens,length($seqS));
        }
        $seqN=$1;
        push(@SeqNs,$seqN);
        $seqS="";
    }
    else{
        $seqS.=$_;
    }
}
$RefSeq{$seqN}=$seqS;
push(@SeqLens,length($seqS));
close $fh_ref;
&VCFConvert($mergefile,$VCFfile);

sub VCFConvert{
    my ($infile,$outfile)=@_;
    open(my $fh_outfile,">$outfile");
    {
my $VCF_Head = << "VCF_HEADER";
##fileformat=VCFv4.1
##source=breakdancer
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=>
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakpoint">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INV,Description="Inverstion">
##ALT=<ID=BND,Description="Breakends">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">
VCF_HEADER
    print $fh_outfile $VCF_Head;
    }

    for(my $i=0;$i<@SeqNs;$i++){
        print $fh_outfile "\#\#contig=\<ID=$SeqNs[$i],length=$SeqLens[$i]\>\n";
    }
    print $fh_outfile "\#\#reference=$reference\n";

    print $fh_outfile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$Sample_Name\n";

    open(my $fh_infile,"$infile");
    my $SV_count=0;
    while(<$fh_infile>){
        chomp();
        my @lines=split(/\t/,$_);
        if($lines[4] ne "ITX" && $lines[4] ne "CTX"){
            my $length=$lines[3]-$lines[1];
            if($lines[4] eq "DEL"){
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4];SVLEN=\-$length\tDP\t$lines[5]\n";
            }
            elsif($lines[4] eq "INV"){
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4];SVLEN=$length\tDP\t$lines[5]\n";
            }
            else{
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4]\tDP\t$lines[5]\n";
            }
        }
        else{
            $SV_count++;
            my $RefA=getRef($lines[0],$lines[1]);
            my $RefB=getRef($lines[2],$lines[3]);
            my $mateidA=$SV_count;
            my $mateidB=$SV_count+1;
            print $fh_outfile "$lines[0]\t$lines[1]\tbnd_$mateidA\t$RefA\t$RefA\]$lines[2]\:$lines[3]\]\t.\tPASS\tSVTYPE=BND;MATEID=bnd_$mateidB\tDP\t$lines[5]\n";
            print $fh_outfile "$lines[2]\t$lines[3]\tbnd_$mateidB\t$RefB\t$RefB\]$lines[0]\:$lines[1]\]\t.\tPASS\tSVTYPE=BND;MATEID=bnd_$mateidA\tDP\t$lines[5]\n";
            $SV_count++;
        }
    }
    close $fh_infile;
    close $fh_outfile;
}

unlink $Tempfile1;
unlink $Tempfile2;

exit(0);