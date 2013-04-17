#!/usr/local/bin/perl

=head1 NAME

    embl_to_fasta.pl

=head1 SYNOPSIS
 
    embl_to_fasta.pl input_embl output_fasta outputdir 
        where input_embl is the input embl file,
              output_fasta is the output fasta file,
              outputdir is the output directory for writing output files. 

=head1 DESCRIPTION

    This script takes an input embl file (<input_embl>), and converts it to 
    fasta format, and writes the output file (<output_fasta>) in directory
    <outputdir>.

=head1 VERSION
  
    Perl script last edited 3-Jan-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script embl_to_fasta.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 3-Jan-13.
# Last edited 3-Jan-2013.
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 3)
{
    print "Usage of embl_to_fasta.pl\n\n";
    print "perl embl_to_fasta.pl <input_embl> <output_fasta> <outputdir>\n";
    print "where <input_embl> is the input embl file,\n";
    print "      <output_fasta> is the output fasta file,\n";
    print "      <outputdir> is the output directory for writing output files\n";
    print "For example, >perl embl_to_fasta.pl Pk_strainH_chr01.embl Pk_strainH_chr01.fa\n";
    print "/lustre/scratch108/parasites/alc/RNA-SeQC/Pknowlesi\n";
    exit;
}

# FIND THE PATH TO THE INPUT EMBL FILE:                     

my $input_embl             = $ARGV[0];

# FIND THE PATH TO THE OUTPUT FASTA FILE:

my $output_fasta           = $ARGV[1];

# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      

my $outputdir              = $ARGV[2];

#------------------------------------------------------------------#

# TEST SUBROUTINES: 

my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_to_output($outputdir);
&test_convert_embl_to_fasta($outputdir);
&test_print_error;

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($outputdir,$input_embl,$output_fasta);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_embl          = $_[1]; # THE INPUT EMBL FILE  
   my $output_fasta        = $_[2]; # THE OUTPUT FASTA FILE 
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 

   # READ IN THE INPUT EMBL FILE, AND MAKE THE OUTPUT FASTA FILE:
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
}

#------------------------------------------------------------------#

# READ IN THE INPUT EMBL FILE, AND MAKE THE OUTPUT FASTA FILE:

sub convert_embl_to_fasta
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN
   my $input_embl          = $_[1]; # INPUT EMBL FILE
   my $output_fasta        = $_[2]; # OUTPUT FASTA FILE
   my $line;                        # 
   my @temp;                        # 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $name                = "none";# NAME OF THE SEQUENCE
   my $expected_length     = "none";# EXPECTED LENGTH OF THE SEQUENCE  
   my $found_seq_start     = 0;     # SAYS WHETHER THE SEQUENCE WAS FOUND 
   my $seq                 = "";    # SEQUENCE
   my $i;                           # 
   my $found_file_end      = 0;     # SAYS WHETHER WE FOUND THE END OF THE FILE
   my $length;                      # LENGTH OF SEQUENCE $seq 

   # READ IN THE INPUT EMBL FILE:
   open(EMBL,"$input_embl") || die "ERROR: convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   while(<EMBL>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\s+/,$line);
      if (substr($line,0,2) eq 'ID')
      {
         if ($name ne 'none' || $expected_length ne 'none') 
         {
            $errormsg      = "ERROR: convert_embl_to_fasta: name $name expected_length $expected_length\n";
            $errorcode     = 1; # ERRORCODE=1 
            return($errorcode,$errormsg); 
         }
         $name             = $temp[1];
         $expected_length  = $temp[$#temp-1];
      } 
      elsif (substr($line,0,2) eq 'SQ')
      {
         if ($found_seq_start == 1)
         {
            $errormsg      = "ERROR: convert_embl_to_fasta: found_seq_start $found_seq_start\n";
            $errorcode     = 3; # ERRORCODE=3 
            return($errorcode,$errormsg);
         }
         $found_seq_start  = 1;
      }
      elsif ($found_seq_start == 1 && substr($line,0,2) ne '//')
      {
         for ($i = 1; $i <= $#temp-1; $i++) { $seq = $seq.$temp[$i];}    
      }
      elsif ($found_seq_start == 1 && substr($line,0,2) eq '//')
      {
         $found_file_end   = 1;
      }
   }
   close(EMBL);
   if ($name eq 'none' || $expected_length eq 'none')
   {
      $errormsg            = "ERROR: convert_embl_to_fasta: name $name expected_length $expected_length\n";
      $errorcode           = 2; # ERRORCODE=2 
      return($errorcode,$errormsg);
   }
   if ($found_seq_start == 0)
   {
      $errormsg            = "ERROR: convert_embl_to_fasta: found_seq_start $found_seq_start\n";
      $errorcode           = 4; # ERRORCODE=4 
      return($errorcode,$errormsg);
   }
   if ($found_file_end == 0)
   {
      $errormsg            = "ERROR: convert_embl_to_fasta: found_file_end $found_file_end\n";
      $errorcode           = 5; # ERRORCODE=5 
      return($errorcode,$errormsg);
   }
   $length                 = length($seq);
   if ($length != $expected_length)
   {
      $errormsg            = "ERROR: convert_embl_to_fasta: length $length expected_length $expected_length\n";
      $errorcode           = 7; # ERRORCODE=7 
      return($errorcode,$errormsg);
   }

   # WRITE OUT THE SEQUENCE IN THE OUTPUT FASTA FILE:
   ($errorcode,$errormsg)  = &print_to_output($output_fasta,$seq,$name);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }     
}

#------------------------------------------------------------------#

# TEST &convert_embl_to_fasta

sub test_convert_embl_to_fasta
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $input_embl;                  # NAME OF EMBL FILE
   my $output_fasta;                # NAME OF OUTPUT FASTA FILE
   my $errorcode           = 0;     # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR 
   my $expected_output_fasta;       # FILE CONTAINING THE EXPECTED CONTENTS OF $output_fasta 
   my $differences;                 # DIFFERENES BETWEEN $output_fasta AND $expected_output_fasta
   my $length_differences;          # LENGTH OF $differences
   my $line;                        # 

   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   system "wget --quiet http://www.maths.tcd.ie/~avrillee/tests/embl_to_fasta_1.embl -O $input_embl";
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test1\n"; exit;}
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_fasta") || die "ERROR: test_convert_embl_to_fasta: cannot open $expected_output_fasta\n";
   print EXPECTED ">SC10H5\n";
   print EXPECTED "gatcagtagacccagcgacagcagggcggggcccagcaggccggccgtggcgtagagcgc\n";
   print EXPECTED "gaggacggcgaccggcgtggccaccgacaggatggctgcggcgacgcggacgacaccgga\n";
   print EXPECTED "gtgtgccagggcccaccacacgccgatggccgcgagcgcgagtcccgcgctgccgaacag\n";
   print EXPECTED "ggcccacagcacactgcgcagaccggcggccacgagtggcgccaggacggtgcccagcag\n";
   print EXPECTED "gagcagcagggtgacgtgggcgcgcgctgcactgtggccgccccgtccgcccgacgcgcg\n";
   print EXPECTED "cggctcgtcatctcgcggtcccaccaccggtcggccccattactcgtcctcaaccctgtg\n";
   print EXPECTED "gcgactgacgttccccggacaggtcgtaccgattgccgccacgccccaccacgcacaggg\n";
   print EXPECTED "cccagacgacgaagcctgacatggtgatcatgacgacggaccacaccgggtagtacggca\n";
   print EXPECTED "gcgagaggaagttggcgatgatcaccagcccggcgatggcgaccccggtgacacgtgccc\n";
   print EXPECTED "acatcgccgttttgagcagcccggcgctgacgaccatggcgagcgcgccgagcgcgagat\n";
   print EXPECTED "ggatccacccccacccggtgagatcgaactggaaaacgtagttgggcgtggtgacgaaga\n";
   print EXPECTED "cgtcgtcctcggcgatggccatgatgccccggaagaggctgagcagcccggcgaggaaga\n";
   print EXPECTED "gcatcaccgccgcgaaggcggtaaggcccgtcgcccattcctgcctcgcggtgtgtgccg\n";
   print EXPECTED "ggtggtgggtatgtgacgtggtcatctcggacctcgtttcgtggaatgcggatgcttcag\n";
   print EXPECTED "cgagcggaggcgccggtgcccgccgcgcccgtgtgccctgccgggccgtgaccggacagg\n";
   print EXPECTED "accaattccttcgccttgcggaactcctcgtccgtgatggcaccccggtctcggatctcg\n";
   print EXPECTED "gagagccgggccagctcgtcgacgctgctggacccgccgcccacggtcttcctgatgtag\n";
   print EXPECTED "gcgtcgaactcctcctgctgagcccgtgcccgcgttgtctcccggctgcccatgttcttg\n";
   print EXPECTED "ccgcgagcgatcacgtagacgaaaacgcccaggaagggcaggaggatgcagaacaccaac\n";
   print EXPECTED "cagccggccttcgcccagccactcagtccgtcgtcccggaagatgtcggtgacgacgcgg\n";
   print EXPECTED "aagagcaggacgaaccacatgatccacaggaagatcatcagcatcgtccagaaggcaccc\n";
   print EXPECTED "agcagtgggtagtcgtacgccaggtaggtctgtgcactcatgtccgtcctccgtcctccg\n";
   print EXPECTED "gggcgcggcccggcggccctcgttccgtactgacatcagggtggtcacgggtcccaccgg\n";
   print EXPECTED "tcggcatcacccggcacgggtgagtggggcgccgaggccgtcgtggtcaggcccgggaca\n";
   print EXPECTED "ccggtgtgaccctggtggaaggacgcgtcccgtggggcacgcaccgccggccgagggcga\n";
   print EXPECTED "ccaccgcctcggtcagtccgagcaggcccagccacaggccgagaagtcgggtcagggcac\n";
   print EXPECTED "gggccgactcggcgggcagcgcgaggacgacgattccggcgacgtcgacggccagcgggt\n";
   print EXPECTED "tgcgcaggcccagcactccggccggggcgcccggcaccagcgtggcgagggccgatgcca\n";
   print EXPECTED "tgagccaggtccaggaacccccaagcctggcgaggacgtgcgccggatcgctcaatgctc\n";
   print EXPECTED "cggtgaccgccccgcccgacccgtctcccttgtcggcaggttccgccgcatcacgcggaa\n";
   print EXPECTED "cggagatggctcccctgtggatcgggcggccgctgcggggccgcccggttggtcggtcgg\n";
   print EXPECTED "tgagcgccggactcccccttcagctcttccagggtcggggtcgacaccgaggtcctggat\n";
   print EXPECTED "cacccgtcaggggtgatccgggcatgccgtcgtggcggtgaggtgggatacgggaacgat\n";
   print EXPECTED "cggcccacgggggaccggacgagacgaagagacgtgagatgagcgatacgaactcgggcg\n";
   print EXPECTED "gcgggcgccaggccgcttccggaccggccccacgtggccgactccctttccgccggcgcg\n";
   print EXPECTED "tggccctggtcgctgtcgcacgtcccctgatcgtcacggtcggtctcgtcaccgcctact\n";
   print EXPECTED "acctgcttcccctggacgagagactcagcgccggcaccctggtgtcgctggtgtgcggac\n";
   print EXPECTED "tgctcgcagtccttctggtgttctgctgggaggtgcgggccatcacgcgctccccgcatc\n";
   print EXPECTED "cgcgtctgagagcgatcgagggcctggccgccacgctggtgctgttcctggtcctcttcg\n";
   print EXPECTED "ccggctcctactacctgctgggtcgctccgcgcccggctccttcagcgagccgctgaaca\n";
   print EXPECTED "ggacggacgcgctgtacttcactctgaccacgttcgccaccgtcggcttcggggacatca\n";
   print EXPECTED "ccgcacgctccgagaccgggcggatcctcacgatggcgcagatgacgggagggctactgc\n";
   print EXPECTED "tcgtcggagtcgccgcccgggtgctggcgagcgcagtgcaggcggggctgcaccgacagg\n";
   print EXPECTED "gccggggaccggcggcatcgccacgctccggtgctgcggaggagccggaggccggaccat\n";
   print EXPECTED "gaccgtacccggtggcttcaccgcctccctgccgccggccgagcgagccgcgtacggcag\n";
   print EXPECTED "gaaggcccgtaaaagggcctcacgttcgtgccacggctggtacgagccggggcagcggcg\n";
   print EXPECTED "gcctgaccccgtcgacctgctggagcgccagtccggcgagcgtgtcccggcactcgtgcc\n";
   print EXPECTED "catccgctacggtcgcatgctggagtcgccgttccgcttctaccgcggtgcggcagcgat\n";
   print EXPECTED "catggcggcggacctggcacccctgcccagcagcggactccaggtgcaattgtgcgggga\n";
   print EXPECTED "cgcgcacccgttgaacttccggctcctggcctcaccggagcgccggctggtcttcgacat\n";
   print EXPECTED "caacgacttcgacgagacgctgcccggccccttcgagtgggacgtcaaacggctggcggc\n";
   print EXPECTED "cggattcgtgatcgcggcccggtcgaacggcttctcgtccaaggaacagaaccgcaccgt\n";
   print EXPECTED "tcgggcctgtgtgcgggcctaccgggagcgcatgagggagttcgccgtcatgccgaccct\n";
   print EXPECTED "ggacatctggtacgcccaggacgacgccgaccacgtacggcaactgctggctacggaggc\n";
   print EXPECTED "cagaggagaagctgagcagcggctcagggacgcggctgcgaaggcccgcacacgcaccca\n";
   print EXPECTED "catgagggcgttcgcgaagctcacccgcgtcacggccgagggccggcgcatcacccccga\n";
   print EXPECTED "cccgccgctgatcaccccactcggcgatctgctcaccgacccggccgaagccggccggga\n";
   print EXPECTED "ggaggaactgcggtccgtcgtgaacggctacgcacggtccctgccgcccgagcgccggca\n";
   print EXPECTED "cctgctgcgtcactaccggcttgtggacatggcgcgcaaggtggtcggcgtcggcagtgt\n";
   print EXPECTED "cggcacccgctgctgggtactgcttctgctcggcagggacgacgacgatcctctgctgct\n";
   print EXPECTED "ccaggccaaggaagcctcggaatcggtgctggcggcccacacgggcggcgaacgctacga\n";
   print EXPECTED "ccatcagggccgcagggtcgtggccggccagcgtctgatccagaccaccggtgacatctt\n";
   print EXPECTED "tctcggctgggcgcgcgtcaccggcttcgacggaaaggcccgggacttctacgtgcgtca\n";
   print EXPECTED "actgtgggactggaagggcgtcgcgcggccggaaaccatggggcccgacctgctctccct\n";
   print EXPECTED "cttcgcccggctgtgcggtgcctgcctggcgagggcccacgcccgttccggtgaccccgt\n";
   print EXPECTED "cgcgctcgccgcgtacctgggcggcagcgaccgcttcgacggcgcgctcaccgagttcgc\n";
   print EXPECTED "ccagtcctacgccgatcagaatgaacgcgaccacgaagctctgctggcggcctgccgctc\n";
   print EXPECTED "cggcagggtcacggccgcccgtttgtgaggccgacccgggaacggccggcgggctggcac\n";
   print EXPECTED "acaccgccgccggtcggcgtcattccggaagctgccgcatctccaggacgcgcaggccca\n";
   print EXPECTED "gcgactggcagcgggtgagcaacccgtacagatgggcctcgtcgatcaccgtgccgaaca\n";
   print EXPECTED "gcacggtctggccggacatgacgacgtgctccagctccgggaacgcgttggccagcgtcc\n";
   print EXPECTED "gtgacaggtgtccctcgacgcggatctcgtagcgcacgagcggtcctttcaccgtaggag\n";
   print EXPECTED "ctcgggacaccgcccggggctccgggtcggacggtgctcttggtgacgagcctgcgcctc\n";
   print EXPECTED "gtcgccctccggtgccctcacccagcacaggtgactccaaccgcagtgtcagtgcctttc\n";
   print EXPECTED "agtgcgtcactgtgatcttgacgacgacgatcaccaggccgagcagtacgttgaccgtcg\n";
   print EXPECTED "cggtgacggccaccagtcgtcgcgaggcgcccgcgcggtgcgccgcggcgacggaccagc\n";
   print EXPECTED "ccacctgaccggcgacggcgacggacagcgccagccacagggtgcccgggacgtccagcc\n";
   print EXPECTED "ccagtacggggctgacggcgatggccgcggccggaggcacggcggccttgacgatcggcc\n";
   print EXPECTED "actcctcgcggcacacacgcagaatcacccgccggtccggagtgtgccgcgcgagacgcg\n"; 
   print EXPECTED "ctccgaacagttcggcgtggacgtgagcgatccagaacaccaagctggtgagcaacagca\n";
   print EXPECTED "gaagaaccagttcggcgcgggggaacgagcccagggtgccggcgccgatcacgacggagg\n"; 
   print EXPECTED "ctgcgagcat\n";
   close(EXPECTED);                                
   $differences            = "";
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";
   system "rm -f $expected_output_fasta";

   # TEST FOR ERRORCODE=1 (ID. APPEARS TWICE IN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 1) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test2\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";

   # TEST FOR ERRORCODE=2 (ID. DOES NOT APPEAR IN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "SQ   Sequence 4870 BP; 769 A; 1717 C; 1693 G; 691 T; 0 other;\n";
   print EMBL "     gatcagtaga cccagcgaca gcagggcggg gcccagcagg ccggccgtgg cgtagagcgc        60\n";
   print EMBL "     gaggacggcg accggcgtgg ccaccgacag gatggctgcg gcgacgcgga cgacaccgga       120\n";
   print EMBL "     gtgtgccagg gcccaccaca cgccgatggc cgcgagcgcg agtcccgcgc tgccgaacag       180\n";
   print EMBL "     ggcccacagc acactgcgca gaccggcggc cacgagtggc gccaggacgg tgcccagcag       240\n";
   print EMBL "     gagcagcagg gtgacgtggg cgcgcgctgc actgtggccg ccccgtccgc ccgacgcgcg       300\n";
   print EMBL "     cggctcgtca tctcgcggtc ccaccaccgg tcggccccat tactcgtcct caaccctgtg       360\n";
   print EMBL "//\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 2) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test3\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";
 
   # TEST FOR ERRORCODE=3 (TWO SEQUENCES IN AN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   print EMBL "SQ   Sequence 4870 BP; 769 A; 1717 C; 1693 G; 691 T; 0 other;\n";
   print EMBL "     gatcagtaga cccagcgaca gcagggcggg gcccagcagg ccggccgtgg cgtagagcgc        60\n";
   print EMBL "     gaggacggcg accggcgtgg ccaccgacag gatggctgcg gcgacgcgga cgacaccgga       120\n";
   print EMBL "     gtgtgccagg gcccaccaca cgccgatggc cgcgagcgcg agtcccgcgc tgccgaacag       180\n";
   print EMBL "     ggcccacagc acactgcgca gaccggcggc cacgagtggc gccaggacgg tgcccagcag       240\n";
   print EMBL "     gagcagcagg gtgacgtggg cgcgcgctgc actgtggccg ccccgtccgc ccgacgcgcg       300\n";
   print EMBL "     cggctcgtca tctcgcggtc ccaccaccgg tcggccccat tactcgtcct caaccctgtg       360\n";
   print EMBL "SQ   Sequence 4870 BP; 769 A; 1717 C; 1693 G; 691 T; 0 other;\n";
   print EMBL "     gatcagtaga cccagcgaca gcagggcggg gcccagcagg ccggccgtgg cgtagagcgc        60\n";
   print EMBL "     gaggacggcg accggcgtgg ccaccgacag gatggctgcg gcgacgcgga cgacaccgga       120\n";
   print EMBL "     gtgtgccagg gcccaccaca cgccgatggc cgcgagcgcg agtcccgcgc tgccgaacag       180\n";
   print EMBL "//\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 3) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test4\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";

   # TEST FOR ERRORCODE=4 (NO SEQUENCE IN AN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 4) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test5\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";

   # TEST FOR ERRORCODE=5 (NO END-OF-FILE SYMBOL IN AN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   print EMBL "SQ   Sequence 4870 BP; 769 A; 1717 C; 1693 G; 691 T; 0 other;\n";
   print EMBL "     gatcagtaga cccagcgaca gcagggcggg gcccagcagg ccggccgtgg cgtagagcgc        60\n";
   print EMBL "     gaggacggcg accggcgtgg ccaccgacag gatggctgcg gcgacgcgga cgacaccgga       120\n";
   print EMBL "     gtgtgccagg gcccaccaca cgccgatggc cgcgagcgcg agtcccgcgc tgccgaacag       180\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 5) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test6\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta";

   # TEST FOR ERRORCODE=7 (WRONG SEQUENCE LENGTH IN AN EMBL FILE):
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_fasta: cannot open input_embl $input_embl\n";
   print EMBL "ID   SC10H5 standard; DNA; PRO; 4870 BP.\n";
   print EMBL "SQ   Sequence 4870 BP; 769 A; 1717 C; 1693 G; 691 T; 0 other;\n";
   print EMBL "     gatcagtaga cccagcgaca gcagggcggg gcccagcagg ccggccgtgg cgtagagcgc        60\n";
   print EMBL "     gaggacggcg accggcgtgg ccaccgacag gatggctgcg gcgacgcgga cgacaccgga       120\n";
   print EMBL "     gtgtgccagg gcccaccaca cgccgatggc cgcgagcgcg agtcccgcgc tgccgaacag       180\n";
   print EMBL "//\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_fasta($outputdir,$input_embl,$output_fasta);
   if ($errorcode != 7) { print STDERR "ERROR: test_convert_embl_to_fasta: failed test7\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_fasta"; 
}

#------------------------------------------------------------------#

# WRITE OUT THE SEQUENCE IN THE OUTPUT FASTA FILE:

sub print_to_output
{
   my $output_fasta        = $_[0]; # OUTPUT FASTA FILE
   my $seq                 = $_[1]; # SEQUENCE
   my $name                = $_[2]; # NAME OF THE SEQUENCE
   my $length;                      # LENGTH OF THE SEQUENCE
   my $offset;                      # 
   my $a_line;                      # A LINE OF THE SEQUENCE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR

   # MAKE AN OUTPUT FASTA FILE:
   open(OUTPUT,">$output_fasta") || die "ERROR: print_to_output: cannot open output_fasta $output_fasta\n";

   print OUTPUT ">$name\n";
   $length                 = length($seq);
   $offset                 = 0;
   while ($offset < $length)
   {
      $a_line              = substr($seq,$offset,60);
      print OUTPUT "$a_line\n";
      $offset              = $offset + 60;
   }

   # CLOSE THE OUTPUT FILE:
   close(OUTPUT);  

   return($errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &print_to_output

sub test_print_to_output
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $output_fasta;                # OUTPUT FILE
   my $seq;                         # SEQUENCE
   my $errorcode;                   # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR 
   my $expected_output_fasta;       # FILE CONTAINING EXPECTED CONTENTS OF $output_fasta
   my $differences;                 # DIFFERENCES BETWEEN $output_fasta AND $expected_output_fasta
   my $line;                        # 
   my $length_differences;          # LENGTH OF $differences

   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   $seq                    = "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG";
   ($errorcode,$errormsg)  = &print_to_output($output_fasta,$seq,"seq1");
   if ($errorcode != 0) { print STDERR "ERROR: test_print_to_output: failed test1\n"; exit;}
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED_OUTPUT,">$expected_output_fasta") || die "ERROR: test_print_to_output: cannot open $expected_output_fasta\n";
   print EXPECTED_OUTPUT ">seq1\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGG\n";
   close(EXPECTED_OUTPUT); 
   $differences            = "";
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_print_to_output: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $output_fasta";
   system "rm -f $expected_output_fasta"; 

}

#------------------------------------------------------------------#

# SUBROUTINE TO MAKE A FILE NAME FOR A TEMPORARY FILE:

sub make_filename
{
   my $outputdir             = $_[0]; # OUTPUT DIRECTORY TO WRITE TEMPORARY FILE NAME TO
   my $found_name            = 0;     # SAYS WHETHER WE HAVE FOUND A FILE NAME YET
   my $filename              = "none";# NEW FILE NAME TO USE 
   my $errorcode             = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg              = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $poss_filename;                 # POSSIBLE FILE NAME TO USE
   my $random_number;                 # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
 
   while($found_name == 0)
   {
      $random_number         = rand();
      $poss_filename         = $outputdir."/tmp".$random_number;
      if (!(-e $poss_filename))
      {
         $filename           = $poss_filename;
         $found_name         = 1;
      } 
   } 
   if ($found_name == 0 || $filename eq 'none')
   {
      $errormsg              = "ERROR: make_filename: found_name $found_name filename $filename\n";
      $errorcode             = 6; # ERRORCODE=6 
      return($filename,$errorcode,$errormsg);
   }

   return($filename,$errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &print_error

sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR

   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 

   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}

#------------------------------------------------------------------#

# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT

   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }

   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }

   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 

   return($errormsg,$errorcode);
}

#------------------------------------------------------------------#

