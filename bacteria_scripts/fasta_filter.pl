#!/usr/bin/perl -w

# Selects only given samples from fasta file

use strict;
use warnings;

my $sample_file = $ARGV[0];
my $fasta_file = $ARGV[1];

open(SAMPLES, $sample_file) || die("Could not open $sample_file: $!\n");
my %sample_list;
while (my $line_in = <SAMPLES>)
{
   chomp $line_in;
   $sample_list{$line_in} = 1;
}

close SAMPLES;

open(FASTA, $fasta_file) || die("Could not open $fasta_file: $!\n");
my $line_in = <FASTA>;
while (!eof(FASTA))
{
   chomp $line_in;

   $line_in =~ m/^>(.+)$/;
   if (!defined($sample_list{$1}))
   {
      while($line_in = <FASTA>)
      {
         if ($line_in =~ m/^>(.+)$/)
         {
            last;
         }
      }
   }
   else
   {
      print $line_in "\n";
      while($line_in = <FASTA>)
      {
         if ($line_in =~ m/^>(.+)$/)
         {
            last;
         }
         else
         {
            print $line_in;
         }
      }
   }
}
close FASTA;

exit(0);


