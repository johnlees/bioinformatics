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

      if ($serotype_list eq "")
      {
         $lane =~ m/^(\d+)_(\d+)(#|_)(\d+)/;
         my $srst2_lane = "$1_$2_$4";
         $serotype_list = `sed '1d' /lustre/scratch108/bacteria/jl11/mlst_serotype/serotype/srst2/results/$srst2_lane\__$lane.95_capsule_sequences_srst2.scores | cut -f 1`;

         chomp $serotype_list;
         my @serotype_lines = split("\n", $serotype_list);
         foreach my $serotype (@serotype_lines)
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
      }
      else
      {
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
      }

      print "\t$format_sero";
   }
   print "\n";
}

exit(0);

