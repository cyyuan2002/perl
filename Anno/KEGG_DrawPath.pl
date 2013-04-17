#!/usr/bin/perl
use strict;
use SOAP::Lite;
use LWP::Simple;
use Getopt::Long; 

my $version="1.0, Alvin Chen, 2011-09-16";

use vars qw($parFileNames $parEnhpercent $parFilterMode $parPvalue $parBgFile $parOutDir $parLogFile);
GetOptions("f=s" =>\$parFileNames,"d=s"=>\$parOutDir,"m:i"=>\$parFilterMode,"b:s"=>\$parBgFile,"v:i"=>\$parPvalue,"l:s"=>\$parLogFile,"p:i" =>\$parEnhpercent);

$parFilterMode=0 if($parFilterMode eq "");
$parEnhpercent=30 if($parEnhpercent eq "");
$parPvalue=0.05 if($parPvalue eq "");
$parLogFile="$parOutDir.log" if($parLogFile eq "");

&Usage if($parFileNames eq "" || $parOutDir eq "");
&Usage if($parFilterMode==2 && $parBgFile eq "");

my $iserror;

my @FileNs=split(/,/,$parFileNames);
foreach my $filename(@FileNs){
  &Filecheck($filename,0);
}
&Filecheck($parBgFile,0) if($parFilterMode==2);

exit(0) if($iserror==1);

my $parWSDL='http://soap.genome.jp/KEGG.wsdl';
my $Soap_ser=SOAP::Lite->service($parWSDL);

my %pathinfors;
my $genecount;



foreach my $filename(@FileNs){
  open(my $fh_file,$filename);
  while(<$fh_file>){
    chomp();
    my @lines=split(/\t/,$_);
    my $pathname=shift(@lines);
    $genecount+=scalar(@lines);
    if(!exists($pathinfors{$pathname})){
      $pathinfors{$pathname}=\@lines;
    }
    else{
      my @genes=@{$pathinfors{$pathname}};
      push(@genes,@lines);
      $pathinfors{$pathname}=\@genes;
    }
  }
  close $fh_file;
}

if(!(-e $parOutDir)){
  mkdir ($parOutDir);
}


if ($parFilterMode<2) {
  open(my $fh_logfile,">$parLogFile");
  foreach my $pathname (keys %pathinfors) {
    my @gnames=@{$pathinfors{$pathname}};
    my @pathgenes=@{$Soap_ser->get_genes_by_pathway($pathname)};
    my $pathnum=scalar(@pathgenes);
    my $gnum=scalar(@gnames);
    #print "$pathname\t",join(",",@gnames),"\t$pathnum\t$gnum\n";
    my $pathpercent;
    if($parFilterMode==1){
      $pathpercent=$gnum*100/$pathnum;
      next if($pathpercent<$parEnhpercent);
    }
    my @fgs;
    my @bgs;
    for (my $i=0;$i<@gnames;$i++) {
      push(@fgs,"#ff0000");
      push(@bgs,"#ffff00");
    }
    my $obj_url=$Soap_ser->color_pathway_by_objects($pathname,\@gnames,\@fgs,\@bgs);
    my @url=split("/",$obj_url);
    my $filename=pop(@url);
    my $pic=get $obj_url;
    $filename="$parOutDir\/$filename";
    open (my $fh_pic,">$filename") || print STDERR "Error: Can't Create File $filename\n";
    print $fh_pic $pic;
    close $fh_pic;
    print $fh_logfile "$pathname\t$gnum\tpathpercent\t",join(",",@gnames),"\n";
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
  open(my $fh_logfile,">$parLogFile");
  foreach my $pathname(keys %pathinfors){
    my @gnames=@{$pathinfors{$pathname}};
    my $gnum=scalar(@gnames);
    if(exists($seqratio{$pathname})){
      my @ratio=split(/:/,$seqratio{$pathname});
      my $pvalue=&fishertest($gnum,$genecount,$ratio[0],$ratio[1]);
      if($pvalue<=$parPvalue){
        my @fgs;
        my @bgs;
        for(my $i=0;$i<@gnames;$i++){
          push(@fgs,"#ff0000");
          push(@bgs,"#ffff00");
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
        print $fh_logfile "$pathname\t$gnum:$genecount\t$seqratio{$pathname}\t$pvalue\t",join(",",@gnames),"\n";
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
		-f     Files of GeneIds (example:fileA,fileB,...)
		-d     Output Directory
		-m     Mode of Filter (0:off, 1:percentage, 2:enrichment)
		-b     Background of the species (must be given when -m 2)
                -v     Pvalue cutoff of fisher-exact test (default:0.05, used when -m 2)
		-p     Percentage cutoff of the Filter (used when -m 1)
                -l     LogFile (default: OutDirectory.log)

    Usage

	exit(0);
};
