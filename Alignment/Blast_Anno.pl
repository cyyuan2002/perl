use strict;
use Getopt::Long;

my $version="1.0 Alvin Chen 2011-04-2";
my %opts;
##Read parameters
GetOptions(\%opts,"i=s","o=s","c:i","help");
if((!defined $opts{i})||(!defined $opts{o})){
    &Usage();
}

my $parInputFile=$opts{i};
my $parOutputFile=$opts{o};
my $parCoverage=(defined $opts{c}) ? $opts{c} : 60;
my $parOutFilter="$parOutputFile.fout";

my $iserror=0;
my @fileInput; #store all the psl information
my %jumplines; ##the line numbers contains redundant infors

&Filecheck($parInputFile,0);

if($iserror==1){
    exit(1);
}

open(fileIN,$parInputFile);
@fileInput=<fileIN>;
close fileIN;

my $fileOut;
open($fileOut,">$parOutputFile");
open(my $filefileter,">$parOutFilter");

my %checkedids;
my %currentinfo;
my $genecount=0;

for(my $i=0;$i<@fileInput;$i++){
    ##read current result information
    my $currline=$fileInput[$i];
    $currline=~s/\n//g;
    %currentinfo=%{&formatblat($currline)};

    #next if($currentinfo{'qsize'}<$parLength);
    #next if($currentinfo{'identity'}<$parIdentity);
    if($currentinfo{'coverage'}<$parCoverage){
	my $qN=$currentinfo{'qname'};
	print $filefileter "$currline\n";
	next;
    }
    next if(exists($checkedids{$currentinfo{'qname'}}));

    $genecount++;
    #my @currgene=@currinfo;

    for(my $j=$i+1;$j<@fileInput;$j++){
	my $newline=$fileInput[$j];
	$newline=~s/\n//g;
	my %newinfo=%{&formatblat($newline)};
	next if ($newinfo{'coverage'}<$parCoverage);
	next if(exists($checkedids{$newinfo{'qname'}}));
	my $isoverlap=&checkoverlap(\%currentinfo,\%newinfo);

	if ($isoverlap==0){ ##no overlap move to next currline
	    my @exons=@{$currentinfo{'exons'}};
	    my $exon;
	    for(my $i=0;$i<@exons;$i++){
		$exon=$exon.$exons[$i]->{'ts'}.":".$exons[$i]->{'te'}.",";
	    }

	    print $fileOut "$currentinfo{'chr'}\t$currentinfo{'chr_S'}\t$currentinfo{'chr_E'}\t$currentinfo{'strand'}\t$currentinfo{'exonnum'}\t$exon\t$currentinfo{'qname'}\tLocus_$genecount\n";
	    last;
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
print $fileOut "$currentinfo{'chr'}\t$currentinfo{'chr_S'}\t$currentinfo{'chr_E'}\t$currentinfo{'strand'}\t$currentinfo{'exonnum'}\t$exon\t$currentinfo{'qname'}\tLocus_$genecount\n";
close $fileOut;
close $filefileter;
exit(0);


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
    $currentinfo{'strand'}=$currentinfo{'strand'}.",".$sinfo->{'strand'};
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
    if($tinfo->{'qsize'}==0){
	my $qN=$tinfo->{'qname'};
	print stderr "Error\t$qN\n";
    }
    my $qmp=$lengthmatch/$tinfo->{'qsize'};
    my $smp=$lengthmatch/$sinfo->{'qsize'};

    if($qmp > 0.3|| $smp > 0.3){
	return 1;
    }
    else{
	return 0;
    }
}

sub formatblat{
    my $blatinfo=shift;
    my @lineinfo=split(/\t/,$blatinfo);
    my @starts=split(/,/,$lineinfo[9]);
    my @exons;
    my %blat;
    for(my $i=0;$i<@starts;$i++){
	my %exoninfo;
	my @SE=split(/:/,$starts[$i]);
	$exoninfo{'ts'}=$SE[0];
	$exoninfo{'te'}=$SE[1];
	push(@exons,\%exoninfo);
    }

    $blat{'chr'}=$lineinfo[3];
    $blat{'chr_S'}=$lineinfo[5];
    $blat{'chr_E'}=$lineinfo[6];
    $blat{'strand'}=$lineinfo[4];
    $blat{'exonnum'}=$lineinfo[7];
    $blat{'qsize'}=$lineinfo[1];
    $blat{'qname'}=$lineinfo[0];
    if($lineinfo[4] eq "-"){
	@exons=reverse(@exons);
    }
    $blat{'exons'}=\@exons;
    $blat{'coverage'}=$lineinfo[2];
    return (\%blat);
}


sub checkoverlap{
    my ($infoA,$infoB)=@_;

    return (0) if($infoA->{'chr'} ne $infoB->{'chr'});  #different chrom
    return (0) if($infoA->{'chr_E'} <= $infoB->{'chr_S'}); # no overlap
    return (1) if($infoA->{'strand'} ne $infoB->{'strand'}); # overlap but different direction
    return (2); #overlap
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
		-c     Coverage of blast hits, default 60

    Usage

	exit(0);
};
