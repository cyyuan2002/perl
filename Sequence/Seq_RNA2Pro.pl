#!/usr/bin/perl
# Translate DNA into protein

use strict;



my $filename=shift;
my $seqName;
my $sequence;
open(in,"$filename");
while(<in>){
	if(/^>(\S+)/){
		$seqName=$1;
	}
	elsif(/(\w+)/){
		$sequence.=$1;
	}
}

my $protein1 = '';
my $protein2 = '';
my $protein3 = '';
my $codon;




# Translate each three-base codon into an amino acid, and append to a protein 
for(my $i=0; $i < (length($sequence) - 2) ; $i += 3) {
    $codon = substr($sequence,$i,3);
    $protein1 .= codon2aa($codon);
}
for(my $i=1; $i < (length($sequence) - 2) ; $i += 3) {
    $codon = substr($sequence,$i,3);
    $protein2 .= codon2aa($codon);
}
for(my $i=2; $i < (length($sequence) - 2) ; $i += 3) {
    $codon = substr($sequence,$i,3);
    $protein3 .= codon2aa($codon);
}

print "$seqName Translate Result\n";
print "Form_1\n";
my $length = length $protein1;   ## bug 1 2004-12-2 14:53
my $modulus = $length % 50;
$protein1 =~ s/(.{50})/$1\n/ig;
#       print  "modulus == $modulus\n";
if ($modulus ==0) {
	print "$protein1";
}
else {
	print "$protein1\n";
}
print "Form_2\n";
$length = length $protein2;   ## bug 1 2004-12-2 14:53
$modulus = $length % 50;
$protein2 =~ s/(.{50})/$1\n/ig;
#       print  "modulus == $modulus\n";
if ($modulus ==0) {
	print "$protein2";
}
else {
	print "$protein2\n";
}
print "From_3\n";
$length = length $protein3;   ## bug 1 2004-12-2 14:53
$modulus = $length % 50;
$protein3 =~ s/(.{50})/$1\n/ig;
#       print  "modulus == $modulus\n";
if ($modulus ==0) {
	print "$protein3";
}
else {
	print "$protein3\n";
}

exit;

sub codon2aa {
    my($codon) = @_;
 
       if ( $codon =~ /GC./i)        { return 'A' }    # Alanine
    elsif ( $codon =~ /TG[TC]/i)     { return 'C' }    # Cysteine
    elsif ( $codon =~ /GA[TC]/i)     { return 'D' }    # Aspartic Acid
    elsif ( $codon =~ /GA[AG]/i)     { return 'E' }    # Glutamic Acid
    elsif ( $codon =~ /TT[TC]/i)     { return 'F' }    # Phenylalanine
    elsif ( $codon =~ /GG./i)        { return 'G' }    # Glycine
    elsif ( $codon =~ /CA[TC]/i)     { return 'H' }    # Histidine
    elsif ( $codon =~ /AT[TCA]/i)    { return 'I' }    # Isoleucine
    elsif ( $codon =~ /AA[AG]/i)     { return 'K' }    # Lysine
    elsif ( $codon =~ /TT[AG]|CT./i) { return 'L' }    # Leucine
    elsif ( $codon =~ /ATG/i)        { return 'M' }    # Methionine
    elsif ( $codon =~ /AA[TC]/i)     { return 'N' }    # Asparagine
    elsif ( $codon =~ /CC./i)        { return 'P' }    # Proline
    elsif ( $codon =~ /CA[AG]/i)     { return 'Q' }    # Glutamine
    elsif ( $codon =~ /CG.|AG[AG]/i) { return 'R' }    # Arginine
    elsif ( $codon =~ /TC.|AG[TC]/i) { return 'S' }    # Serine
    elsif ( $codon =~ /AC./i)        { return 'T' }    # Threonine
    elsif ( $codon =~ /GT./i)        { return 'V' }    # Valine
    elsif ( $codon =~ /TGG/i)        { return 'W' }    # Tryptophan
    elsif ( $codon =~ /TA[TC]/i)     { return 'Y' }    # Tyrosine
    elsif ( $codon =~ /TA[AG]|TGA/i) { return '_' }    # Stop
    else {
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}