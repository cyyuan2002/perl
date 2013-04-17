#!/usr/bin/perl
use strict;
use SOAP::Lite;
use LWP::Simple;
my ($parFilein,$parOutDir)=@ARGV;
if(@ARGV<2){
  print "Usage:$0 InputFile OutputDir\n";
  exit(0);
}
my $parWSDL='http://soap.genome.jp/KEGG.wsdl';
my $Soap_ser=SOAP::Lite->service($parWSDL);

my %paths;
open(my $fh_filein,$parFilein) || die "Can't open file $parFilein\n";
while (<$fh_filein>) {
  chomp();
  my @lines=split(/\t/,$_);
  my $gname=$lines[0];
  my @kegg_info=split(/\s+/,$Soap_ser->bconv("ncbi-geneid:$gname"));
  my $pathway=$Soap_ser->get_pathways_by_genes([$kegg_info[1]]);
  if ($pathway) {
    foreach my $hit (@{$pathway}) {
      if(exists($paths{$hit})){
        my $refgenames=$paths{$hit}->{gnames};
        my @genames=@{$refgenames};
        push(@genames,$lines[0]);
        $paths{$hit}->{gnames}=\@genames;
      }
      else{
        my @genames;
        push(@genames,$lines[0]);
        $paths{$hit}->{gnames}=\@genames;
      }
    }
  }
}
close $fh_filein;

mkdir ($parOutDir);

foreach my $key(keys %paths){
  my %tempinfo=%{$paths{$key}};
  my @gnames=@{$tempinfo{gnames}};
  next if(scalar(@gnames)<3);
  #print join(",",@gnames),"\n";
  my @fgs;
  my @bgs;
  for(my $i=0;$i<@gnames;$i++){
    push (@fgs,"#ff0000");
    push (@bgs,"#ffff00");
  }
  my $obj_url=$Soap_ser->color_pathway_by_objects($key,\@gnames,\@fgs,\@bgs);
  my @url=split("/",$obj_url);
  my $filename=pop(@url);
  my $pic=get $obj_url;
  $filename="$parOutDir\/$filename";
  open (my $fh_pic,">$filename") || print STDERR "Error: Can't Create File $filename\n";
  print $fh_pic $pic;
  close $fh_pic;
  print "$obj_url\n";
}


sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

