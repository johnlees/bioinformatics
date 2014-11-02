#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 ) . "/../assembly_scripts/";

use assembly_common;

my $help_message = <<HELP;
Usage: ./mutate_fasta.pl -f <fasta_file> -r <mutation_rate> -o <output_file> (-s <seed>) > mutations.txt

Creates SNPs in a multifasta using JC69 rate matrix
Options:
   -f --fasta    Multi-fasta file of sequence to create SNPs in
   -r --rate     Per-base mutation rate
   -s --seed     Optional seed for random number generator
   -o --output   Output file name

   -h --help     This help message

Prints mutations tab separated to stdout and logs to stderr
HELP

my $total_mutations;

#
# Main
#
my ($fasta_file, $mutation_rate, $output_file, $seed, $help);
GetOptions ("fasta|f=s"  => \$fasta_file,
            "rate|r=f" => \$mutation_rate,
            "output|o=s"  => \$output_file,
            "seed|s=i"  => \$seed,
            "help|h"     => \$help
		   ) or die($help_message);

# Parse input
if (defined($help))
{
   print $help_message;
}
elsif (!defined($fasta_file) || !defined($mutation_rate) || !defined($output_file))
{
   print STDERR "One of the required options was not defined\n";
   print $help_message;
}
elsif (!-e $fasta_file)
{
   print STDERR "The file $fasta_file does not exist\n";
}
else
{
   # Options ok, pass sequences to mutation routine
   my $multi_fasta = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta') || die("Failed to open $fasta_file: $!\n");
   my $fasta_out = Bio::SeqIO->new(-file => ">$output_file", -format => 'fasta') || die("Couldn't write to $output_file: $!\n");

   # Loop through sequences
   while (my $sequence = $multi_fasta->next_seq())
   {
      my $seq_string = $sequence->seq();
      my ($new_seq, $mutations);

      if (defined($seed))
      {
         ($new_seq, $mutations) = assembly_common::create_snps($seq_string, $mutation_rate, $seed);
      }
      else
      {
         ($new_seq, $mutations) = assembly_common::create_snps($seq_string, $mutation_rate, $seed);
      }

      $sequence->seq($new_seq);
      $fasta_out->write_seq($sequence);

      # Add mutations to list
      print join("\t", @{$mutations}) . "\n";
      $total_mutations += scalar(@$mutations);
   }

   # Give total number of mutations actually made
   print STDERR "$total_mutations mutations produced\n";

   # Done
   print STDERR "\nDone\n";
}

exit(0);

