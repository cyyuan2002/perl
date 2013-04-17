#!/usr/bin/perl
use strict;

my ($parInfile,$parFasfile,$parTempdir,$parOutfile)=@ARGV;
if(@ARGV<3){
    print "Usage:$0 Input_file Fasta_file Temp_dir [Output_file]\n";
    exit(1);
}
my $iserror=0;
&Filecheck($parInfile,0);
&Filecheck($parFasfile,0);
&Filecheck($parTempdir,0);

my %seqs;
if($parOutfile eq ""){
    $parOutfile="$parInfile.contigs";
}
my $parLogfile="$parOutfile.seqlog";

if($iserror==1){
    exit(1);
}

open(my $fh_fasfile,$parFasfile);
my $seqName;
my $seq;
while(<$fh_fasfile>){
    chomp();
    if(/^>(\S+)/){
	if($seq ne ""){
	    $seqs{$seqName}=$seq;
	}
	$seqName=$1;
	$seq="";
    }
    else{
	$seq=$seq.$_;
    }
}
$seqs{$seqName}=$seq;
close $fh_fasfile;

open(my $fh_outfile,">$parOutfile");
open(my $fh_outlog,">$parLogfile");
open(my $fh_outdir,">$parOutfile.dir");
open(my $fh_infile,$parInfile);
while(<$fh_infile>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[2]!=2){
	my $contigN=$lines[1];
	my $seqids=$lines[0];
	my @fasseqs=&getseq($seqids);
	&dealassembly($lines[1],@fasseqs);
    }
}

close $fh_infile;
close $fh_outlog;
close $fh_outfile;
exit(0);

sub dealassembly{
    my $printresult=0;
    my ($locusid,@fasseqs)=@_;
#    if($locusid eq "Locus_128"){
#	$printresult=1;
#    }
    if(scalar(@fasseqs)<2){
	my $seqN=$fasseqs[0]->{"id"};
	my $seq=$fasseqs[0]->{"seq"};
	print $fh_outlog "\#$locusid\n>$seqN\n$seq\n";
	#print "DDDD>$locusid\_0\n$seq\n";
	print $fh_outfile ">$locusid\_0\n$seq\n";
	print $fh_outdir "$locusid\_0\t$seqN\t+\n";
    }
    else{
	my $fastaseqs;
	for(my $i=0;$i<@fasseqs;$i++){
	    my $seqN=$fasseqs[$i]->{"id"};
	    my $seq=$fasseqs[$i]->{"seq"};
	    $fastaseqs.=">$seqN\n$seq\n";
	}
	my $tempfile="$parTempdir\/temp.fas";
	open(my $fh_tempfile,">$tempfile");
	print $fh_tempfile $fastaseqs;
	my $result=`cap3 $parTempdir\/temp.fas`;
	#print "AAAA\n$result" if($printresult==1);
	my @contigstrand=&detectstrand($result); # get strand info from assembly information

	my $contigfile="$parTempdir\/temp.fas.cap.contigs";
	open(my $confile,$contigfile) || die "Can't open file $contigfile";
	my @contigs=<$confile>;
	my $seqnum=scalar(@contigstrand);
	#print "seqnum:$seqnum\n" if($printresult==1);
	my $conseq=&dealseqid($locusid,@contigs);
	#my $contigstr=join("",@contigs);
	#my $seqnum=$contigstr=~s/>/>/g;
	#print "$locusid\tSeqNum1:\t";
	if($seqnum>1){
	    my @fasseqs2=&checkfas(@fasseqs);
	    my $isNs=shift(@fasseqs2);
	    if($isNs==0){
		print $fh_outlog "\#$locusid\n$fastaseqs\n";
		print $fh_outfile "$conseq";
		for(my $i=0;$i<@contigstrand;$i++){
		    print $fh_outdir "$locusid\_$i\t$contigstrand[$i]\n";
		}
	    }
	    else{
		my $fastaseqs2;
		for(my $i=0;$i<@fasseqs2;$i++){
		    my $seqN=$fasseqs2[$i]->{"id"};
		    my $seq=$fasseqs2[$i]->{"seq"};
		    $fastaseqs2.=">$seqN\n$seq\n";
		}
		my $tempfile="$parTempdir\/temp.fas";
		open(my $fh_tempfile,">$tempfile");
		print $fh_tempfile $fastaseqs2;
		my $result2=`cap3 $parTempdir\/temp.fas`;
		#print "$fastaseqs2" if($printresult==1);
		#print "$result2" if($printresult==1);
		my @contigstrand2=&detectstrand($result2);
		my $contigfile2="$parTempdir\/temp.fas.cap.contigs";
		open(my $confile2,$contigfile2);
		my @contigs2=<$confile2>;
		my $conseq2=&dealseqid($locusid,@contigs2);
		my $seqnum2=scalar(@contigstrand2);
		#print "SeqNum2:$seqnum2\n";
		if($seqnum<=$seqnum2){
		    print $fh_outlog "\#$locusid\n$fastaseqs\n";
		    print $fh_outfile "$conseq";
		    for(my $i=0;$i<@contigstrand;$i++){
			print $fh_outdir "$locusid\_$i\t$contigstrand[$i]\n";
		    }
		}
		else{
		    print $fh_outlog "\#$locusid\n$fastaseqs2\n";
		    print $fh_outfile "$conseq2";
		    for(my $i=0;$i<@contigstrand2;$i++){
			print $fh_outdir "$locusid\_$i\t$contigstrand2[$i]\n";
		    }
		}
	    }
	}
	elsif($seqnum==1){
	    my $conseq=&dealseqid($locusid,@contigs);
	    print $fh_outlog "\#$locusid\n$fastaseqs\n";
	    #print "CCCC$conseq";
	    print $fh_outfile "$conseq";
	    print $fh_outdir "$locusid\_0\t$contigstrand[0]\n";
	}
	else{ ##seqnum=0 no assembly
	    print $fh_outlog "\#$locusid\n$fastaseqs\n";
	    for(my $i=0;$i<@fasseqs;$i++){
		my $seqN=$fasseqs[$i]->{"id"};
		my $seq=$fasseqs[$i]->{"seq"};
		print $fh_outfile ">$locusid\_$i\n$seq\n";
		print $fh_outdir "$locusid\_$i\t$seqN\t+\n";

	    }
	}
	my $command="$parTempdir\/temp.*";
	`rm $command`;
    }
}

sub detectstrand{
    my $assinfo=shift;
    my @assinfos=split(/\n/,$assinfo);
    my @dirinfo;
    my $seqids;
    my $dirs;
    my $ismatch=0;
    #print "$assinfo\n";
    for(my $i=0;$i<@assinfos;$i++){
	if($assinfos[$i]=~/\*+\s(\S+\s\d+)\s\*+/){
	    if($seqids eq ""){
		$ismatch=1;
		next;
	    }
	    else{
		my $dirinfo="$seqids\t$dirs";
		push(@dirinfo,$dirinfo);
	    }
	    $seqids="";
	    $dirs="";
	}
	else{
	    if($ismatch==1){
		if($assinfos[$i]=~/(\S+)([\+\-])/){
		    if($seqids eq ""){
			$seqids=$1;
			$dirs=$2;
		    }
		    else{
			$seqids.=",$1";
			$dirs.=",$2";
		    }
		}
		elsif($assinfos[$i]=~/\s+(\S+)([\+\-])/){
		    if($seqids eq ""){
			$seqids=$1;
			$dirs=$2;
		    }
		    else{
			$seqids.=",$1";
			$dirs.=",$2";
		    }
		}
		elsif($assinfos[$i]=~/DETAILED DISPLAY OF CONTIGS/){
		    last;
		}
	    }
	}

    }
    if($seqids ne ""){
	my $dirinfo="$seqids\t$dirs";
	push(@dirinfo,$dirinfo);
    }
    return @dirinfo;
}

sub dealseqid{
    my @conseq=@_;
    my $locusid=shift(@conseq);
    my $count=0;
    my $fasseqs;
    for(my $i=0;$i<@conseq;$i++){
	if($conseq[$i]=~/^>\S+/){
	    $fasseqs.=">$locusid\_$count\n";
	    $count++;
	}
	else{
	    $fasseqs.="$conseq[$i]";
	}
    }

    return $fasseqs;
}

sub getseq{
    my $seqids=shift;
    my @seqid=split(/,/,$seqids);
    my @seqfastas;
    for(my $i=0;$i<@seqid;$i++){
	if(exists($seqs{$seqid[$i]})){
	    #$seq=&checkfas($seqid[$i],$seqs{$seqid[$i]});
	    my %seqinfo=("id"=>$seqid[$i],"seq"=>$seqs{$seqid[$i]});
	    push(@seqfastas,\%seqinfo);
	}
	else{
	    print "Error:AAA\n";
	    exit(1);
	}
    }
    return (@seqfastas);
}

sub checkfas{
    my @seqinfors=@_;
    my @newseqs;
    my $isNs=0;
    for(my $i=0;$i<@seqinfors;$i++){
	my $seqN=$seqinfors[$i]->{"id"};
	my $seq=$seqinfors[$i]->{"seq"};
	my @seqment;
	my $count;
	while($seq=~/(N{5,})/){
	    $isNs=1;
	    my $end=$+[0];
	    my $length=length($1);
	    my $start=$end-$length;
	    my $seq1=substr($seq,0,$start);
	    $seq=substr($seq,$end);
	    push (@seqment,$seq1);
	}
	push(@seqment,$seq);
	if(@seqment<2){
	    my %seqinfo=("id"=>$seqN,"seq"=>$seq);
	    push(@newseqs,\%seqinfo);
	}
	else{
	    #print "AAAAA\n";
	    for(my $j=0;$j<@seqment;$j++){
		my %seqinfo=("id"=>"$seqN\_$j","seq"=>$seqment[$j]);
		#print ">$seqN\_$j\n$seqment[$j]\n";
		push(@newseqs,\%seqinfo);
	    }
	    #print "BBBB\n";
	}
    }
    #print "IsNs\t$isNs\n";
    return ($isNs,@newseqs);
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
