#!/usr/bin/perl -w

use strict;
use warnings;

my @columns = ("hsdM,1.1", "hsdM,1.2", "1.1,2.1", "1.1,2.2", "1.1,2.3", "1.2,2.1", "1.2,2.2", "1.2,2.3");

open(SAMPLES, "../strep_pairs_lanes.txt") || die("Could not open lanes\n");

while (my $sample = <SAMPLES>)
{
   chomp $sample;
   $sample =~ /^(\d+_\d+)[#_](\d+)$/;

   my $sample_name = "$1_$2";

   if (!open(IVR, "$sample_name/ivr_mapped_allele.txt"))
   {
      print STDERR ("Could not open $sample_name output\n");
   }
   else
   {
      my %reads_mapped;
      my $ivr_header = <IVR>;
      while (my $ivr_line = <IVR>)
      {
         chomp $ivr_line;
         my ($upstream, $D39, $R6, $reads) = split("\t", $ivr_line);

         $reads_mapped{$upstream}{$D39} = $reads;
      }

      print $sample;

      my (%max, %max_name);
      foreach my $column (@columns)
      {
         my ($upstream, $map) = split(',', $column);

         if (!defined($reads_mapped{$upstream}{$map}))
         {
            $reads_mapped{$upstream}{$map} = 0;
         }
         elsif (!defined($max{$upstream}) || $reads_mapped{$upstream}{$map} > $max{$upstream})
         {
            $max{$upstream} = $reads_mapped{$upstream}{$map};
            $max_name{$upstream} = $map;
         }

         print "\t$reads_mapped{$upstream}{$map}";
      }

      if (defined($max_name{"hsdM"}) && defined($max_name{$max_name{"hsdM"}}))
      {
         print "\t" . $max_name{'hsdM'} . "-" . $max_name{$max_name{'hsdM'}} ."\n";
      }
      else
      {
         print "\t\n";
      }

      close IVR;
   }

}

close SAMPLES;

exit(0);

