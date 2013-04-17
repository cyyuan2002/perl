#!/usr/bin/perl
#===============================================================================
#
#         FILE:KEGG_Regulated_Drawpath.pl
#
#        USAGE:
#
#  DESCRIPTION:This program is used to draw KEGG pathway with up/down regulated gene information
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.0
#      CREATED: 2012-07-21
#     REVISION:
#===============================================================================

use strict;
use SOAP::Lite;
use LWP::Simple;
use Getopt::Long;

my $version="1.0, Alvin Chen, 2012-07-21";

use vars qw($parFileName $parSpecies $parEnhnum $parFilterMode $parPvalue $parBgFile $parOutDir $parLogFile);
GetOptions("f=s" =>\$parFileName,"s=s" => \$parSpecies,"d=s"=>\$parOutDir,"m:i"=>\$parFilterMode,"b:s"=>\$parBgFile,"v:i"=>\$parPvalue,"l:s"=>\$parLogFile,"n:i" =>\$parEnhnum);

$parFilterMode ||=0;
$parEnhnum ||=2;
$parPvalue ||= 0.05;
$parLogFile ||="$parOutDir.log";
$parSpecies ||= "map";

&Usage if($parFileName eq "" || $parOutDir eq "");
&Usage if($parFilterMode==2 && $parBgFile eq "");

my $iserror;

my @FileNs=split(/,/,$parFileName);
foreach my $filename(@FileNs){
  &Filecheck($filename,0);
}
&Filecheck($parBgFile,0) if($parFilterMode==2);

exit(0) if($iserror==1);

my $parWSDL='http://soap.genome.jp/KEGG.wsdl';
my $Soap_ser=SOAP::Lite->service($parWSDL);

my %pathinfos;
my %upgenes;
my %downgenes;
my %genelst;
my $genecount;

#----File Example-------
#CNAG_00393      CNA03810        2.4.1.18        U
#CNAG_06000      CNM00090        2.7.1.- U
#CNAG_05932      CNF00430        5.2.1.8 D
#CNAG_04189      CNI03270        1.3.5.1 D
#---------------------


open(my $fh_filein,$parFileName) || die "Can't open file $parFileName\n";
while (<$fh_filein>) {
  chomp();
  $genecount++;
  my @lines=split(/\t/,$_);
  my $ncbiID="$parSpecies:$lines[1]";
  $genelst{$ncbiID}="$lines[0]|$lines[2]";
  if($lines[3] eq "U"){
    $upgenes{$ncbiID}=1;
  }
  else{
    $downgenes{$ncbiID}=1;
  }

  my $pathway=$Soap_ser->get_pathways_by_enzymes([$lines[2]]);
  if ($pathway) {
    foreach my $hit (@{$pathway}) {
        $hit=~s/map/$parSpecies/g;
        if(exists($pathinfos{$hit})){
            my $refgenames=$pathinfos{$hit};
            my @genames=@{$refgenames};
            push(@genames,$ncbiID);
            $pathinfos{$hit}=\@genames;
        }
        else{
            my @genames;
            push(@genames,$ncbiID);
            $pathinfos{$hit}=\@genames;
        }
    }
  }
}
close $fh_filein;

open(my $fh_logfile,">$parLogFile");

#foreach my $key(keys %pathinfos){
#  my @gnames=@{$pathinfos{$key}};
#  print $fh_logfile "$key\t";
#  foreach my $key(@gnames){
#    print $fh_logfile "$genelst{$key},";
#  }
#  print $fh_logfile "\n";
#}

if(!(-e $parOutDir)){
  mkdir ($parOutDir);
}

if ($parFilterMode<2) {
  open(my $fh_logfile,">$parLogFile");
  foreach my $pathname (keys %pathinfos) {
    my @gnames=@{$pathinfos{$pathname}};
    #my @pathgenes=@{$Soap_ser->get_genes_by_pathway($pathname)};
    #my $pathnum=scalar(@pathgenes);
    my $gnum=scalar(@gnames);
    #print "$pathname\t",join(",",@gnames),"\t$pathnum\t$gnum\n";
    if($parFilterMode==1){
      next if($gnum<$parEnhnum);
    }
    my @fgs;
    my @bgs;
    for (my $i=0;$i<@gnames;$i++) {
        if(exists($upgenes{$gnames[$i]})){
            push(@fgs,"yellow");
            push(@bgs,"red");
        }
        else{
            push(@fgs,"green");
            push(@bgs,"blue");
        }
    }
    my $obj_url=$Soap_ser->color_pathway_by_objects($pathname,\@gnames,\@fgs,\@bgs);
    my @url=split("/",$obj_url);
    my $filename=pop(@url);
    my $pic=get $obj_url;
    $filename="$parOutDir\/$filename";
    open (my $fh_pic,">$filename") || print STDERR "Error: Can't Create File $filename\n";
    print $fh_pic $pic;
    close $fh_pic;
    print $fh_logfile "$pathname\t$gnum\t";
    foreach my $key(@gnames){
        print $fh_logfile "$genelst{$key},";
    }
    print $fh_logfile "\n";
  }
  close $fh_logfile;
}
else{
  open (my $fh_bgfile,$parBgFile);
  my %seqratio;
  while(<$fh_bgfile>){
    chomp();
    my @lines=split(/\t/,$_);
    $seqratio{$lines[0]}=$lines[1];
  }
  close $fh_bgfile;
  foreach my $pathname(keys %pathinfos){
    my @gnames=@{$pathinfos{$pathname}};
    my $gnum=scalar(@gnames);
    if(exists($seqratio{$pathname})){
      my @ratio=split(/:/,$seqratio{$pathname});
      my $pvalue=&fishertest($gnum,$genecount,$ratio[0],$ratio[1]);
      if($pvalue<=$parPvalue){
        my @fgs;
        my @bgs;
        for(my $i=0;$i<@gnames;$i++){
            if(exists($upgenes{$gnames[$i]})){
                push(@fgs,"yellow");
                push(@bgs,"red");
            }
            else{
                push(@fgs,"green");
                push(@bgs,"blue");
            }
        }
        my $obj_url=$Soap_ser->color_pathway_by_objects($pathname,\@gnames,\@fgs,\@bgs);
        my @url=split("/",$obj_url);
        my $filename=pop(@url);
        my $pic=get $obj_url;
        $filename="$parOutDir\/$filename";
        open (my $fh_pic,">$filename") || print STDERR "Error: Can't Create File $filename\n";
        print $fh_pic $pic;
        close $fh_pic;
        $pvalue=sprintf("%.4f",$pvalue);
        print $fh_logfile "$pathname\t$gnum:$genecount\t$seqratio{$pathname}\t$pvalue\t";
        foreach my $key(@gnames){
            print $fh_logfile "$genelst{$key},";
        }
        print $fh_logfile "\n";
      }
    }
    else{
      print STDERR "Can't find pathway: $pathname\n";
    }
  }
  close $fh_logfile;
}


sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub fishertest{
  my ($n11,$n12,$n21,$n22)=@_;
  use Text::NSP::Measures::2D::Fisher2::right;
  my $n1p=$n11+$n12;
  my $np1=$n11+$n21;
  my $npp=$n11+$n12+$n21+$n22;
  my $right_value=calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
  return $right_value;
}

sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print STDERR "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print STDERR "Directory $Filename already exists\n";
	    $iserror=1;
	}
    }
}


sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-f     File of Up/Down regulated genes
                -s     KEGG spcices ID (example: mmu,cne)
		-d     Output Directory
		-m     Mode of Filter (0:off, 1:number, 2:enrichment)
		-b     Background of the species (must be given when -m 2)
                -v     Pvalue cutoff of fisher-exact test (used when -m 2, default: 0.05)
		-n     Number cutoff of the Filter (used when -m 1, default: 2)
                -l     LogFile (default: OutDirectory.log)

    Usage

	exit(0);
};
