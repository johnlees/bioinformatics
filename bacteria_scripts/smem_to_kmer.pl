#!/usr/bin/perl -w

use strict;
use warnings;

my $kmer_file = $ARGV[0];
my $smem_file = $ARGV[1];

open(KMER, $kmer_file) || die("Could not open $kmer_file\n");

my @kmers;
while (my $line_in = <KMER>)
{
   my $line_in = <KMER>;
   chomp $line_in;

   push(@kmers, $line_in);
}

close KMER;

open(SMEM, $smem_file) || die("Could not open $smem_file\n");

my $i = -1;
while (my $line_in = <SMEM>)
{
   $i++;
   chomp $line_in;
   my ($junk1, $junk2, $length) = split("\t", $line_in);

   while ($line_in = <SMEM>)
   {
      chomp $line_in;
      if ($line_in eq "//")
      {
         last;
      }
      elsif ($line_in eq "")
      {
         next;
      }
      else
      {
         my ($junk1, $start, $end, @info) = split("\t", $line_in);
         if ($start == 0 && $end == $length)
         {
            print "$kmers[$i] 1\n";
            while ($line_in = <SMEM>)
            {
               chomp $line_in;
               if ($line_in eq "//")
               {
                  last;
               }
            }
            last;
         }
      }
   }

}


exit(0);

