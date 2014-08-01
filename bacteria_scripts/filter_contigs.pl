#!/usr/bin/perl -w

##
# filter_contigs.pl
# Script to remove low coverage and short contigs/scaffolds from de novo
# assemblies
####

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;

# Default cutoffs
my $cov_cutoff = 3;
my $len_cutoff = 300;

my $help_message = <<END;
Usage: filter_contigs.pl -i <input_file> -o <output_file> <options>
Removes low coverage and short contigs from a de novo assembly

	-i, --input    Name of input multi-fasta contig file
	-o, --output   Name of output multi-fasta

	-c, --cov      Coverage cutoff. Default 3x
	-l, --len      Length cutoff. Default 300bp

	-h, --help     Displays this help message
END

# Get input files and parameters
my ($input_file, $output_file, $help);
GetOptions ("input|i=s"  => \$input_file,
            "output|o=s" => \$output_file,
            "cov|c=o"  => \$cov_cutoff,
            "len|l=o"    => \$len_cutoff,
            "help|h"     => \$help
		   ) or die($help_message);

# Check necessary files exist
if (defined($help))
{
   print $help_message;
}
elsif (!defined($input_file) || !-e $input_file)
{
	print ("Input file does not exist!\n\n");
	print $help_message;
}
else
{
   # Do filtering
   my $fasta_in = Bio::SeqIO->new(-file => $input_file,
                                  -format => 'fasta') || die($!);

   my $fasta_out = Bio::SeqIO->new(-file => ">$output_file",
                                  -format => 'fasta') || die($!);

   while (my $contig = $fasta_in->next_seq())
   {
      my ($length, $coverage);
      if ($contig->id =~ /^NODE_\d+_length_(\d+)_cov_(.+)_ID_\d+$/)
      {
         $length = $1;
         $coverage = $2;
      }

      if ($length >= $len_cutoff && $coverage >= $cov_cutoff)
      {
         $fasta_out->write_seq($contig);
      }
   }
}

exit(0);

