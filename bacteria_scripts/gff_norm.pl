#!/usr/bin/perl -w

#
# Normalises gff files
# i.e. turns format of 23FSpn reference into R6
#
# Usage:
# ./gff_norm.pl 23FSpn_noseq.gff > new.gff
#

use strict;
use warnings;

my $input_gff = $ARGV[0];

open(GFF, $input_gff) || die("Couldn't open $input_gff");

while (my $line_in = <GFF>)
{
   chomp $line_in;

   if ($line_in =~ /^#/)
   {
      print $line_in . "\n";
   }
   else
   {
      my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split("\t", $line_in);

      if ($feature eq "CDS")
      {
         my @attrs = split(";", $attribute);

         my ($new_id, $gene);
         foreach my $att_pair (@attrs)
         {
            my ($key, $value) = split("=", $att_pair);

            if ($key eq "ID")
            {
               if ($value =~ /"(.+)"/)
               {
                  $new_id = "\"$1.1\"";
               }
               else
               {
                  $new_id = "\"$value.1\"";
               }
            }
            elsif ($key eq "gene")
            {
               $gene = $value
            }
         }

         print join("\t", $seqname, $source, "gene", $start, $end, ".", $strand, $frame, "ID=$new_id;gene=$gene") . "\n";
      }
      print $line_in . "\n";
   }
}

exit(0);

