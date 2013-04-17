#!/usr/bin/perl
use strict;
use SOAP::Lite;
use LWP::Simple;
my ($parFilein,$parFileout)=@ARGV;
if(@ARGV<2){
  print "Usage:$0 Inputfile Outputfile\n";
  exit(0);
}

my $parWSDL='http://soap.genome.jp/KEGG.wsdl';
my $Soap_ser=SOAP::Lite->service($parWSDL);

my %paths;
open(my $fh_filein,$parFilein) || die "Can't open file $parFilein\n";
my $genecount=0;
while (<$fh_filein>) {
  chomp();
  $genecount++;
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
open (my $fh_fileout,">$parFileout");

foreach my $key(keys %paths){
  my %tempinfo=%{$paths{$key}};
  my @gnames=@{$tempinfo{gnames}};
  my @pathgenes=@{$Soap_ser->get_genes_by_pathway($key)};
  my $pathnum=scalar(@pathgenes);
  my $gnum=scalar(@gnames);
  print $fh_fileout "$key\t$gnum:$genecount\t$pathnum\t",join(",",@gnames),"\t",join(",",@pathgenes),"\n";
}
close $fh_fileout;

exit(1);


sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

