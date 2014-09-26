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

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use assembly_pipeline;

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
   assembly_pipeline::filter_contigs($len_cutoff, $cov_cutoff, $input_file, $output_file);
}

exit(0);

