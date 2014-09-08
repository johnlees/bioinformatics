#!/usr/bin/perl -w

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use Getopt::Long;

use quake_wrapper;

my $help_message = <<HELP;
Usage: ./run_quake.pl -r <reads_file> -k <kmer_size> -t <threads>

Runs quake error correction on paired end reads provided as fastq(.gz) files.
Options:
   -r --reads    Tab delimited file of fastq or fastq.gz file locations.
                 One line per sample. First column sample name, second
                 column forward reads, third column reverse reads.
                 It is best to give absolute paths
   -k --kmer     The kmer size to use with quake. Set to
                 ceil[log(200*G)/log(4)] where G is the genome length
   -t --threads  The number of threads to allow quake to use

   -h --help     This help message

Output will be in cwd/quake
HELP

#
# Main
#
my ($read_file, $kmer_size, $threads, $help);
GetOptions ("kmer|k=i"  => \$kmer_size,
            "threads|t=s" => \$threads,
            "reads|r=s"  => \$read_file,
            "help|h"     => \$help
		   ) or die($help_message);

# Parse input
if (defined($help))
{
   print $help_message;
}
elsif (!defined($read_file) || !defined($kmer_size) || !defined($threads))
{
   print STDERR "One of the options was not defined. All are required\n";
   print $help_message;
}
elsif (!-e $read_file)
{
   print STDERR "The file $read_file does not exist\n";
}
else
{
   # Options ok, run quake

   # Get read locations in a hash
   my($samples, $reads) = quake_wrapper::parse_read_file($read_file);

   # Run quake
   quake_wrapper::quake_error_correct($reads, $kmer_size, $threads);

   # Done
   print "\nDone\n";
}

exit(0);

