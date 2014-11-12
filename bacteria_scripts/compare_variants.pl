#!/usr/bin/perl -w

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use Getopt::Long;

use compare_variants;

my $help_message = <<HELP;
Usage: ./compare_variants.pl --vcf1 <1.vcf.gz> --vcf2 <2.vcf.gz> --ref1 <ref1.fa> --ref2 <ref2.fa>

Compares variant calls in vcfs that have been mapped/placed on different reference
sequences

   Options

   Required
   --vcf1              First vcf file. Can be vcf/bcf/bgzipped
   --vcf2              Second vcf file.
   --ref1              Reference fasta that vcf1 variants have been placed on
   --ref2              Reference fasta that vcf2 variants have been placed on


   Optional
   --strict            Requires ref and alt positions to match as well as
                       sequence surrounding them. May fail for INDELs
   --threshold         Minimum blast score to claim a match. Default 280
   --top-hit           Report only the top blast hit for each variant in vcf1

   -h, --help          Shows more detailed help.

Requires: bcftools, blastn
HELP

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($vcf1, $ref1, $vcf2, $ref2, $strict, $threshold, $top_hit, $help);
GetOptions( "vcf1=s" => \$vcf1,
            "vcf2=s" => \$vcf2,
            "ref1=s" => \$ref1,
            "ref2=s" => \$ref2,
            "strict" => \$strict,
            "top-hit" => \$top_hit,
            "threshold=i" => \$threshold,
            "help|h" => \$help
		   ) or die($help_message);

if (defined($help))
{
   print STDERR "$help_message";
}
elsif (!-e $vcf1 || !-e $vcf2 || !-e $ref1 || !-e $ref2)
{
   print STDERR "ERROR: One of the four required input files does not exist\n";
}
else
{
   my @vcfs = ($vcf1, $vcf2);
   my @refs = ($ref1, $ref2);

   # Set threshold for score
   if (!defined($threshold))
   {
      $threshold = 280;
   }

   # Run blastn comparison of windows
   compare_variants::create_blastn_input(\@vcfs, \@refs, "blast_windows");
   my $blast_scores = compare_variants::blastn_pairwise("blast_windows.1.fa", "blast_windows.2.fa");

   # Print blast results
   foreach my $q_id (sort keys %$blast_scores)
   {
      my @hits;
      my $match = 0;
      my $total_matches = scalar(keys %{$$blast_scores{$q_id}});

      foreach my $s_id (sort keys %{$$blast_scores{$q_id}})
      {
         $match++;
         if ($strict)
         {
            my ($q_chrom, $q_pos, $q_ref, $q_alt) = split(',', $q_id);
            my ($s_chrom, $s_pos, $s_ref, $s_alt) = split(',', $s_id);

            # Skip over if ref and alt do not match (only really works for
            # SNPs)
            if (!(($q_ref eq $s_ref && $q_alt eq $s_alt) || ($q_ref eq $s_alt && $q_alt eq $s_ref)))
            {
               next;
            }
         }

         # Only print scores above provided threshold
         if ($$blast_scores{$q_id}{$s_id} >= $threshold)
         {
            my $output_string = "$q_id\t$s_id\t$$blast_scores{$q_id}{$s_id}\n";

            # Only print top hit if required (store in array until the last
            # match reached)
            if ($top_hit && $total_matches > 1)
            {
               push(@hits, $output_string);

               # Reached last match. Print the one with the highest match score
               if ($match == $total_matches)
               {
                  my $best_hit;
                  my $high_score = 0;

                  foreach my $hit (@hits)
                  {
                     chomp $hit;
                     my ($q, $s, $score) = split("\t", $hit);

                     if ($score > $high_score)
                     {
                        $high_score = $score;
                        $best_hit = $hit;
                     }
                  }
                  print "$best_hit\n";
               }
            }
            # Otherwise print everything
            else
            {
               print $output_string;
            }
         }
      }
   }

   unlink "blast_windows.1.fa", "blast_windows.2.fa";
}

exit(0);

