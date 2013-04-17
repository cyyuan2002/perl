#!/usr/bin/env perl
###############################
# gff2gtf.pl
#
# parses mRNA, exon lines from a gff file and prints gtf lines (for cufflinks)
# 5/2010
# #############################

use lib "/home/SCE/chenyuan/bin/Perl5Lib/lib/perl5";
use Bio::Tools::GFF;

my $parser = new Bio::Tools::GFF->new(-file=> $ARGV[0], -gff_version => 3);

my %hash;
while( my $result = $parser->next_feature )
{
	($id,@junk)= $result->get_tag_values("ID");
	$type = $result->primary_tag();

	if(!$result)
	{
		last;
	}

	$seq_id = $result->seq_id();
	$strand = $result->strand();
	$strand =~ s/-1/-/g;
	$strand =~ s/1/+/g;
	$start = $result->start();
	$end = $result->end();

	if($type eq "mRNA")
	{
		($parent,@junk)= $result->get_tag_values("Parent");
		$hash{$id} = $parent;
	}
	if($type eq "exon")
	{
		#find out transcript (parent) and gene for THIS exon
		($parent,@junk)= $result->get_tag_values("Parent");
		$transcript = $parent;
		$gene = $hash{$transcript};
		print "$seq_id\tFlyBase\t$type\t$start\t$end\t.\t$strand\t.\tgene_id \"$gene\";transcript_id \"$transcript\";\n";
	}
}
