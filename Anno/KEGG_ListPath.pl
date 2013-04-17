#!/usr/bin/perl
use strict;
use SOAP::Lite;
use LWP::Simple;
my $parFilein=shift;

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

foreach my $key(keys %paths){
  my %tempinfo=%{$paths{$key}};
  my @gnames=@{$tempinfo{gnames}};
  print "$key\t";
  foreach my $gN(@gnames){
    print "$gN\t";
  }
  print "\n";
}


sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}
