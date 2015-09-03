#!/usr/bin/perl -w

use strict;
use warnings;

use Text::CSV;

sub extract_cn($)
{
   my ($cn) = @_;

   my $int_cn = 1;
   if ($cn =~ /^CN(\d+)$/)
   {
      $int_cn = $1;
   }

   return $int_cn;
}

sub compare_calls($$)
{
   my ($cn1, $cn2) = @_;

   my $int_cn1 = extract_cn($cn1);
   my $int_cn2 = extract_cn($cn2);

   my $diff = 0;
   if ($int_cn1 != $int_cn2 && ($int_cn1 <= 1 || $int_cn2 <= 1))
   {
      $diff = 1;
   }

   return $diff;
}

# Take CSV file as input
my $in_file = $ARGV[0];

# Set up csv reader
my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();

open my $fh, "<:encoding(utf8)", "$in_file" or die "$in_file: $!";

# Get and set header
my $header_row = $csv->getline ($fh);

my $i = 0;
my (%cols, %diff);
foreach my $sample (@$header_row)
{
   if ($sample =~ /^X(\d+)_(.+)$/)
   {
      $cols{$1}{$2} = $i;
      $diff{$1} = 0;
   }
   $i++;
}

# Read each row
my @cnv;
my $num_diff = 0;
$i = 0;

while ( my $row = $csv->getline( $fh ) )
{
   foreach my $sample (keys %cols)
   {
      if (compare_calls($row->[$cols{$sample}{csf}], $row->[$cols{$sample}{blood}]))
      {
         $cnv[$i]++;

         if (!$diff{$sample})
         {
            $num_diff++;
            $diff{$sample} = 1;
         }
      }
   }
   $i++;
}

print "Number of samples different: $num_diff\n";
print "Number of differences per CNV\n";

for ($i = 0; $i < scalar(@cnv); $i++)
{
   if (defined($cnv[$i]) && $cnv[$i] > 0)
   {
      print join("\t", $i, $cnv[$i]) . "\n";
   }
}

exit(0);

