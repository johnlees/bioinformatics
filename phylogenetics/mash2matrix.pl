#!/usr/bin/perl -w

use strict;
use warnings;

my $dist_file = $ARGV[0];
open(INPUT, $dist_file) || die("Could not open $dist_file\n");

my @mat_out;
while (my $line_in = <INPUT>)
{
   chomp $line_in;
   my ($el1, $el2, $dist, @junk) = split("\t", $line_in);

   $el1 =~ m/SE0+(\d+)_contigs\.fa$/;
   my $row = $1 - 1;

   $el2 =~ m/SE0+(\d+)_contigs\.fa$/;
   my $col = $1 - 1;

   $mat_out[$row][$col] = $dist;
}

for (my $i = 0; $i < scalar(@{$mat_out[0]}); $i++)
{
   print join("\t", @{$mat_out[$i]}) . "\n";
}

exit(0);

