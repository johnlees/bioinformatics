#!/usr/bin/perl -w

use strict;
use warnings;

open(META, "metadata.txt") || die("Couldn't open metadata.txt");

while (my $meta_line = <META>)
{
   chomp $meta_line;

   my ($sample, @lanes) = split("\t", $meta_line);
   print $sample;

   foreach my $lane (@lanes)
   {
      my $serotype_list = `grep '$lane\\b' /lustre/scratch108/bacteria/jl11/mlst_serotype/serotype/srst2/nl_s_pneumo_serotype__compiledResults.txt | cut -f 2-`;
      chomp $serotype_list;

      my $format_sero = "";
      my (@serotypes) = split("\t", $serotype_list);

      foreach my $serotype (@serotypes)
      {
         if ($serotype ne "-")
         {
            if ($format_sero ne "")
            {
               $format_sero .= "/$serotype";
            }
            else
            {
               $format_sero = $serotype;
            }
         }
      }

      print "\t$format_sero";
   }
   print "\n";
}

exit(0);

