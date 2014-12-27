#!/usr/bin/perl -w
#
use strict;
use warnings;

use Getopt::Long;
use Text::CSV;

# Parses qc files from bacterial sequencing
# At the moment set up to extract run, lane and tag of all paired isolates

my $usage = "parse_qc_metadata.pl qc_file.csv > pairs.txt\nparse_qc_metadata.pl qc_file.csv --strain 'strain_regex' > pairs.txt\n";

#
# Main
#

my $metadata_in = $ARGV[0];

# Use a regular expression for strain name compared to sample name if specified
my ($strain_regex_in, $strain_regex);
GetOptions("strain=s" => \$strain_regex_in);
if (defined($strain_regex_in))
{
   $strain_regex = qr/$strain_regex_in/;
}

if (!defined($metadata_in) || !-e ($metadata_in))
{
   # Error/usage message
   print "Input file not found\n";
   print $usage;
}
else
{
   my $csv = Text::CSV->new ( { binary => 1, sep_char => ",", auto_diag => 1, eol => $/ } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();
   open my $metadata_csv, "<:encoding(utf8)", "$metadata_in" or die "$metadata_in: $!";

   my %lanes;
   my @pairs;

   # Throw away header
   my $header = $csv->getline($metadata_csv);

   # Read through file, parse fields
   while (my $line_in = $csv->getline($metadata_csv))
   {
      my ($sample, $sanger_id, $mx_pool, $run, $lane, $tag, $total_yield, $ref_matches, $matching_yield, $yield_size_ratio, $strain_name) = @$line_in;

      # Use strain regex if defined
      if ((!defined($strain_regex)) || ($strain_name =~ $strain_regex))
      {
         # Standardise CSF sample names
         if ($sample =~ /^(\d+)_I$/)
         {
            $sample = $1;
         }

         # Put into hash of run_lane_tag vs sample id
         $lanes{$sample} = "$run" . "_$lane#$tag";

         # Look for any blood isolates
         if ($sample =~ /^(\d+)_II$/)
         {
            push(@pairs, $1);
         }
      }
   }

   close $metadata_csv;

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

