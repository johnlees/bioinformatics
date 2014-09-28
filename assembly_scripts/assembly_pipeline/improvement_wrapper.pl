#!/usr/bin/perl -w

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use assembly_pipeline;

my $threads = 4;
my $cov_cutoff = 3;
my $len_cutoff = 300;

my $forward_reads = $ARGV[0];
my $reverse_reads = $ARGV[1];

my $contigs = $ARGV[2];

my $output_dir = $ARGV[3];

if (!defined($forward_reads) || !defined($reverse_reads) || !defined($output_dir) || !defined($contigs))
{
   print STDERR "This is a wrapper script and should not be run directly";
}
else
{
   assembly_pipeline::run_improvement($contigs, $forward_reads, $reverse_reads, $output_dir);
}

exit(0);

