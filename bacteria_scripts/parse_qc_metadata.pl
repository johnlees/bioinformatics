#!/usr/bin/perl -w
#
use strict;
use warnings;

# Parses qc files from bacterial sequencing
# At the moment set up to extract run, lane and tag of all paired isolates

my $usage = "parse_qc_metadata.pl qc_file.csv > pairs.txt\n";

my $metadata_in = $ARGV[0];

if (!defined($metadata_in) || !-e ($metadata_in))
{
   # Error/usage message
   print "Input file not found\n";
   print $usage;
}
else
{
   open(META, "$metadata_in") || die("Could not open $metadata_in");

   my %lanes;
   my @pairs;

   # Throw away header
   my $header = <META>;

   # Read through file, parse fields
   while (my $line_in = <META>)
   {
      chomp($line_in);

      my ($sample, $sanger_id, $mx_pool, $run, $lane, $tag, $total_yield, $ref_matches, $matching_yield, $yield_size_ratio, $strain_name) = split(",", $line_in);

      # Put into hash of run_lane_tag vs sample id
      $lanes{$sample} = "$run" . "_$lane#$tag";

      # Look for any blood isolates
      if ($sample =~ /^(\d+)_II$/)
      {
         push(@pairs, $1);
      }
   }

   close META;

   # Print out lane, run and tag of each pair to stdout
   print "Sample_ID\tCSF_lane\tblood_lane\n";

   foreach my $paired_sample (@pairs)
   {
      my $csf_run = $lanes{$paired_sample};
      my $blood_run = $lanes{$paired_sample . "_II"};

      if (!defined($csf_run))
      {
         print STDERR "Sample $paired_sample has no csf run recorded\n";
      }
      elsif (!defined($blood_run))
      {
         print STDERR "Sample $paired_sample has no blood run recorded\n";
      }
      else
      {
         print join("\t", $paired_sample, $csf_run, $blood_run) . "\n";
      }
   }
}

exit(0);

