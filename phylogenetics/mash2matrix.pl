#!/usr/bin/perl -w

use strict;
use warnings;

my $dist_file = $ARGV[0];
open(INPUT, $dist_file) || die("Could not open $dist_file\n");

my %row_idx;
my @mat_out;
my $i = 0;
while (my $line_in = <INPUT>)
{
   chomp $line_in;
   my ($el1, $el2, $dist, @junk) = split("\t", $line_in);

   if (!defined($row_idx{$el1}))
   {
      $row_idx{$el1} = $i++;
   }
   if (!defined($row_idx{$el2}))
   {
      $row_idx{$el2} = $i++;
   }

   $mat_out[$row_idx{$el1}][$row_idx{$el2}] = $dist;
}

for (my $i = 0; $i < scalar(@{$mat_out[0]}); $i++)
{
   print join("\t", @{$mat_out[$i]}) . "\n";
}

exit(0);

