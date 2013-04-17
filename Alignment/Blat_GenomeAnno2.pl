use strict;
use Getopt::Long;

#my $version="1.0 Alvin Chen 2011-01-16";
my $version="1.1 Alvin Chen 2011-03-28";
my %opts;
##Read parameters
GetOptions(\%opts,"i=s","o=s","a:i","c:i","l:i","help");
if((!defined $opts{i})||(!defined $opts{o})){
    &Usage();
}

my $parInputFile=$opts{i};
my $parOutputFile=$opts{o};
my $parIdentity=(defined $opts{a}) ? $opts{a} : 90;
my $parCoverage=(defined $opts{c}) ? $opts{c} : 90;
my $parLength=(defined $opts{l}) ? $opts{l} : 150;
my $parBaseDiff=5;


my $iserror=0;
my @fileInput; #store all the psl information
my %jumplines; ##the line numbers contains redundant infors


&Filecheck($parInputFile,0);
&Filecheck($parOutputFile,1);

if($iserror==1){
    exit(1);
}

open(fileIN,$parInputFile);
@fileInput=<fileIN>;
close fileIN;

my $fileOut;
open($fileOut,">$parOutputFile");

my %checkedids;
my %currentinfo;
my $genecount=0;

for(my $i=0;$i<@fileInput;$i++){
    ##read current result information
    my $currline=$fileInput[$i];
    $currline=~s/\n//g;
    %currentinfo=%{&formatblat($currline)};
    next if($currentinfo{'qsize'}<$parLength);
    next if($currentinfo{'identity'}<$parIdentity);
    next if(exists($checkedids{$currentinfo{'qname'}}));

    $genecount++;
    #my @currgene=@currinfo;

    for(my $j=$i+1;$j<@fileInput;$j++){
	my $newline=$fileInput[$j];
	$newline=~s/\n//g;
	my %newinfo=%{&formatblat($newline)};
	next if ($newinfo{'qsize'}<$parLength);
	next if ($newinfo{'identity'}<$parIdentity);
	next if(exists($checkedids{$newinfo{'qname'}}));
	my $isoverlap=&checkoverlap(\%currentinfo,\%newinfo);

	if ($isoverlap==0){ ##no overlap move to next currline
	    my @exons=@{$currentinfo{'exons'}};
	    my $exon;
	    for(my $i=0;$i<@exons;$i++){
		$exon=$exon.$exons[$i]->{'ts'}.":".$exons[$i]->{'te'}.",";
	    }

	    print $fileOut "$currentinfo{'chr'}\t$currentinfo{'chr_length'}\t$currentinfo{'chr_S'}\t$currentinfo{'chr_E'}\t$currentinfo{'strand'}\t$currentinfo{'exonnum'}\t$exon\t$currentinfo{'qname'}\tgene_$genecount\n";
	    last;
	}
	elsif ($isoverlap==1){ ##overlap but different direction
	    next;
	}
	else{ ##overlap same direction
	    my $isoform=&checkisoform(\%currentinfo,\%newinfo);
	    if($isoform==1){
		$checkedids{$newinfo{'qname'}}=1;  #isoform conformed
		&mergeinfo(\%currentinfo,\%newinfo);
	    }
	}
    }

}

my @exons=@{$currentinfo{'exons'}};
my $exon;
for(my $i=0;$i<@exons;$i++){
    $exon=$exon.$exons[$i]->{'ts'}.":".$exons[$i]->{'te'}.",";
}

#print $fileOut ""
print $fileOut "$currentinfo{'chr'}\t$currentinfo{'chr_length'}\t$currentinfo{'chr_S'}\t$currentinfo{'chr_E'}\t$currentinfo{'strand'}\t$currentinfo{'exonnum'}\t$exon\t$currentinfo{'qname'}\tgene_$genecount\n";


sub mergeinfo{
    my ($tinfo,$sinfo)=@_;
    my @exons=(@{$tinfo->{'exons'}},@{$sinfo->{'exons'}});
    @exons=&sortexons(@exons);
    my @newexons;
    my $lastS=$exons[0]->{'ts'};
    my $lastE=$exons[0]->{'te'};
    for(my $i=1;$i<@exons;$i++){
	if($exons[$i]->{'ts'} > $lastE){
	    my %exoninfo;
	    $exoninfo{'ts'}=$lastS;
	    $exoninfo{'te'}=$lastE;
	    push(@newexons,\%exoninfo);
	    $lastS=$exons[$i]->{'ts'};
	    $lastE=$exons[$i]->{'te'};
	}
	elsif($exons[$i]->{'ts'}<=$lastE && $exons[$i]->{'te'}>$lastE){
	    $lastE=$exons[$i]->{'te'};
	}
	elsif($exons[$i]->{'ts'}>=$lastS && $exons[$i]->{'te'} <= $lastE){
	    next;
	}
    }
    my %exoninfo;
    $exoninfo{'ts'}=$lastS;
    $exoninfo{'te'}=$lastE;
    push(@newexons,\%exoninfo);

    $currentinfo{'exons'}=\@newexons;
    my ($chrS,$chrE,$arraylength,$querylength)=&chromregion(@newexons);
    $currentinfo{'chr_S'}=$chrS;
    $currentinfo{'chr_E'}=$chrE;
    $currentinfo{'exonnum'}=$arraylength;
    $currentinfo{'qsize'}=$querylength;
    $currentinfo{'qname'}=$currentinfo{'qname'}.",".$sinfo->{'qname'};
}

sub chromregion{
    my @exons=@_;
    my %exon1=%{$exons[0]};
    my $chrS=$exon1{'ts'};
    my $arraylength=scalar(@exons);
    my %exonE=%{$exons[$arraylength-1]};
    my $chrE=$exonE{'te'};
    my $querylength;
    for(my $i=0;$i<@exons;$i++){
	my $start=${$exons[$i]}{'ts'};
	my $end=${$exons[$i]}{'te'};
	$querylength=$querylength+($end-$start+1);
	if($start<$chrS){
	    $chrS=$start;
	}
	if($end>$chrE){
	    $chrE=$end;
	}
    }
    return ($chrS,$chrE,$arraylength,$querylength);
}


sub sortexons{
    my @exons=@_;
    my @sorted;
    while(scalar(@exons)>0){
	my $minindex=0;
	my $minstart=$exons[0]->{'ts'};
	for(my $i=1;$i<@exons;$i++){
	    my $start=$exons[$i]->{'ts'};
	    if($start<$minstart){
		$minindex=$i;
		$minstart=$start;
	    }
	}
	push(@sorted,$exons[$minindex]);
	splice(@exons,$minindex,1);
    }
    return @sorted;
}


sub checkisoform{  ##check the isoform of different blat hit
    my ($tinfo,$sinfo)=@_;
    my @texons=@{$tinfo->{'exons'}};
    my @sexons=@{$sinfo->{'exons'}};

    my $exonmatch=0;
    my $lengthmatch=0;
    my $exonnum=0;
    my $querylength=0;

    my $texoncount=0;
    my $sexoncount=0;

    while($texoncount<@texons && $sexoncount<@sexons){
	my %texonregion=%{$texons[$texoncount]};
	my %sexonregion=%{$sexons[$sexoncount]};

	#my @texonregion=split(/\:/,$texons[$texoncount]);
	#my @sexonregion=split(/\:/,$sexons[$sexoncount]);
	if($texonregion{'ts'}<=$sexonregion{'ts'} && $texonregion{'te'}>=$sexonregion{'ts'}){
	    if($texonregion{'te'}<=$sexonregion{'te'}){
	        my $overlength=$texonregion{'te'}-$sexonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		$texoncount++;
	    }
	    else{
		my $overlength=$sexonregion{'te'}-$sexonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		$sexoncount++;
	    }
	}
	elsif($sexonregion{'ts'}<=$texonregion{'ts'} && $sexonregion{'te'}>=$texonregion{'ts'}){
	    if($texonregion{'te'}<=$sexonregion{'te'}){
		my $overlength=$texonregion{'te'}-$texonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		$texoncount++;
	    }
	    else{
		my $overlength=$sexonregion{'te'}-$texonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		$sexoncount++;
	    }
	}
	else{
	    if($texonregion{'ts'}>$sexonregion{'te'}){
		$sexoncount++;
	    }
	    elsif($sexonregion{'ts'}>$texonregion{'te'}){
		$texoncount++;
	    }
	}
    }

    my $qmp=$lengthmatch/$tinfo->{'qsize'};
    my $smp=$lengthmatch/$sinfo->{'qsize'};

    if($qmp > 0.5|| $smp > 0.5){
	return 1;
    }
    else{
	return 0;
    }
}



sub formatblat{
    my $blatinfo=shift;
    my @lineinfo=split(/\t/,$blatinfo);
    my @sizes=split(/,/,$lineinfo[18]);
    my @starts=split(/,/,$lineinfo[20]);
    my @exons;
    my %blat;
    for(my $i=0;$i<@sizes;$i++){
	my %exoninfo;
	my $start=$starts[$i]+1;
	my $end=$starts[$i]+$sizes[$i];
	$exoninfo{'ts'}=$start;
	$exoninfo{'te'}=$end;
	push(@exons,\%exoninfo);
    }

    my $hitcover=(($lineinfo[0]+$lineinfo[1]+$lineinfo[2])/$lineinfo[10])*100;
    my $pid=get_pid(@lineinfo);

    $blat{'chr'}=$lineinfo[13];
    $blat{'chr_length'}=$lineinfo[14];
    $blat{'chr_S'}=$lineinfo[15]+1;
    $blat{'chr_E'}=$lineinfo[16];
    $blat{'strand'}=$lineinfo[8];
    $blat{'exonnum'}=$lineinfo[17];
    $blat{'qsize'}=$lineinfo[10];
    $blat{'qname'}=$lineinfo[9];
    $blat{'exons'}=\@exons;
    $blat{'coverage'}=$hitcover;
    $blat{'identity'}=$pid;
    return (\%blat);
}


sub checkoverlap{
    my ($infoA,$infoB)=@_;

    return (0) if($infoA->{'chr'} ne $infoB->{'chr'});  #different chrom
    return (0) if($infoA->{'chr_E'} <= $infoB->{'chr_S'}); # no overlap
    return (1) if($infoA->{'strand'} ne $infoB->{'strand'}); # overlap but different direction
    return (2); #overlap
}


sub get_pid {
    my @line = @_;
    my $pid = (100.0 - (&pslCalcMilliBad(@line) * 0.1));
    return $pid;
}


sub pslCalcMilliBad {
    my @cols = @_;

    # sizeNul depens of dna/Prot
    my $sizeMul=1;

    my $qAliSize = $sizeMul * ($cols[12] - $cols[11]);
    my $tAliSize = $cols[16] - $cols[15];

    # I want the minimum of qAliSize and tAliSize
    my $aliSize;
    $qAliSize < $tAliSize ? $aliSize = $qAliSize : $aliSize = $tAliSize;

    # return 0 is AliSize == 0
    return 0 if ($aliSize <= 0);

    # size diff
    my $sizeDiff = $qAliSize - $tAliSize;
    if ($sizeDiff < 0) {
	    $sizeDiff = 0;
    }

    # insert Factor
    my $insertFactor = $cols[4];
    my $milliBad = (1000 * ($cols[1]*$sizeMul + $insertFactor + &round(3*log( 1 + $sizeDiff)))) / ($sizeMul * ($cols[0] + $cols[2] + $cols[1]));
    return $milliBad;
}

sub round {
    my $number = shift;
    return int($number + .5);
}


sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "File $Filename already exists\n";
	    $iserror=1;
	}
    }
}

sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-i     Input psl file
		-o     Output annonation file
		-a     Identity of blat hits, default 90
		-c     Coverage of blat hits, default 90
		-l     Minimal length of target sequence, default 150

    Usage

	exit(0);
};
