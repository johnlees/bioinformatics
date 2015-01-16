#!perl -w

use strict;
use warnings;

my $pair_file = $ARGV[0];
open (PAIRS, $pair_file) || die("Could not open input $pair_file\n");

my $discordant_samples;
my %discordant_alleles;
my %discordant;

my $wc = `wc -l $pair_file`;
my ($null, $lines, $file) = split(/\s+/, $wc);
my $num_samples = ($lines - 3)/6;

my $header = <PAIRS>;
while (my $line_in = <PAIRS>)
{
   chomp $line_in;

   my ($sample_allele, $L95, $U95) = split(/\s+/, $line_in);

   my ($sample, $allele);
   if ($sample_allele =~ /^(\d+)_([A-F])$/)
   {
      $sample = $1;
      $allele = $2;

      if (($L95 < 0 && $U95 > 0) || ($L95 > 0 && $U95 < 0))
      {
         if (!defined($discordant{$sample}))
         {
            $discordant_samples++;
            $discordant{$sample} = 1;
         }
         $discordant_alleles{$allele}++;
      }
   }
   else
   {
      print STDERR "$sample_allele doesn't look like sample_allele\n";
   }
}

close PAIRS;

print "Allele\tDiscordant pairs\tPercentage\n";

my $discordant_percentage = $discordant_samples/$num_samples;
print "All\t$discordant_samples\t$discordant_percentage" . "\n";

foreach my $alleles (sort keys %discordant_alleles)
{
   my $allele_percentage = $discordant_alleles{$alleles}/$num_samples;
   print "$alleles\t$discordant_alleles{$alleles}\t$allele_percentage\n";
}

exit(0);

