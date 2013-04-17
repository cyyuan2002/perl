#!/usr/bin/perl
use strict;
use SOAP::Lite;
use LWP::Simple;
my ($parFilein,$parSpecies)=@ARGV;

die "Usage:$0 <file_in> <KEGG_species>\n" if(@ARGV<2);

#CNAG_04926      CNM02600        1.1.1.100
#CNAG_03555      CNG00470        1.1.1.101
#CNAG_07749      CNG00270        1.1.1.102
#CNAG_00831      CNA08100        1.1.1.103

open(my $fh_filein,$parFilein) || die "Can't open file $parFilein\n";

my %enzymes;
my %paths;
my $genecount=0;
while (<$fh_filein>) {
  chomp();
  $genecount++;
  my @lines=split(/\t/,$_);
  if(!exists($enzymes{$lines[2]})){
    my @genes;
    push(@genes,$lines[0]);
    $enzymes{$lines[2]}=\@genes;
  }
  else{
    my @genes=@{$enzymes{$lines[2]}};
    push(@genes,$lines[0]);
    $enzymes{$lines[2]}=\@genes;
  }
}
close $fh_filein;

my $parWSDL='http://soap.genome.jp/KEGG.wsdl';
my $Soap_ser=SOAP::Lite->service($parWSDL);

foreach my $key(keys %enzymes){
 my $pathway=$Soap_ser->get_pathways_by_enzymes([$key]);
  if ($pathway) {
    foreach my $hit (@{$pathway}) {
      $hit=~s/map/$parSpecies/g;
      if(exists($paths{$hit})){
        my @genames=@{$paths{$hit}};
        push(@genames,@{$enzymes{$key}});
        $paths{$hit}=\@genames;
      }
      else{
        my @genames;
        push(@genames,@{$enzymes{$key}});
        $paths{$hit}=\@genames;
      }
    }
  }
}

foreach my $key(keys %paths){
  my @gnames=@{$paths{$key}};
  my $gnum=scalar(@gnames);
  print "$key\t$gnum:$genecount\t",join(",",@gnames),"\n";
}


sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}
