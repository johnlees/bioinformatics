#!/usr/bin/perl -w

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use assembly_pipeline;

my $contigs = $ARGV[0];
my $sample = $ARGV[1];
my $genus = $ARGV[2];
my $threads = $ARGV[3];

if (!defined($contigs) || !defined($sample) || !defined($genus))
{
   print STDERR "This is a wrapper script and should not be run directly";
}
else
{
   assembly_pipeline::run_annotation($contigs, $sample, $genus, $threads);
}

exit(0);

