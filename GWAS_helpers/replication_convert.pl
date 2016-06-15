#!/usr/bin/perl -w

use strict;
use warnings;

# Converting Thomas Benfield's genotype data
# Usage: ./replication_convert.pl HumanOmni1-Quad_v1-0_H.csv cph_gwas_key.csv gwas_illumina_cph.csv
#

sub conv_call($$$)
{
   my ($call, $rsid, $mapping) = @_;

   my $formatted = "0\t0"; # Missing call default
   if (defined($$mapping{$rsid}{A}))
   {
      if ($call eq "A_A")
      {
         $formatted = join("\t", $$mapping{$rsid}{A}, $$mapping{$rsid}{A});
      }
      elsif ($call eq "A_B")
      {
         $formatted = join("\t", $$mapping{$rsid}{A}, $$mapping{$rsid}{B});
      }
      elsif ($call eq "B_B")
      {
         $formatted = join("\t", $$mapping{$rsid}{B}, $$mapping{$rsid}{B});
      }
   }
   else
   {
      print STDERR "could not find $rsid in manifest\n";
   }

   return $formatted;
}

my $manifest_file = $ARGV[0];
my $metadata_file = $ARGV[1];
my $genotype_file = $ARGV[2];

my $snp_regex = qr/^\[(.)\/(.)\]$/;

open(METADATA, $metadata_file) || die("Could not open metadata file $metadata_file: $!\n");
my $header = <METADATA>;

my %metadata;
while (my $line_in = <METADATA>)
{
   chomp $line_in;

   my ($sample, $pheno, $sex) = split(";", $line_in);
   $metadata{$sample}{pheno} = $pheno + 1;
   $metadata{$sample}{sex} = $sex + 1;
}

close METADATA;

open(MANIFEST, $manifest_file) || die("Could not open manifest file $manifest_file: $!\n");
# header
for (my $i = 1; $i<=7; $i++)
{
   $header = <MANIFEST>;
}

my %mapping;
while (my $line_in = <MANIFEST>)
{
   chomp $line_in;

   my ($name, $rsid, $strand, $snp, $junk1,$junk2,$junk3,$junk4,$junk5,$chr,$pos,@junks) = split(",", $line_in);
   $mapping{$rsid}{chr} = $chr;
   $mapping{$rsid}{pos} = $pos;

   $snp =~ $snp_regex;
   $mapping{$rsid}{A} = $1;
   $mapping{$rsid}{B} = $2;
}

close MANIFEST;

# rsids are in the header, otherwise file is ped-ish
open(GENOTYPES, $genotype_file) || die("Could not open genotype file $genotype_file: $!\n");
$header = <GENOTYPES>;
chomp $header;
my @rsids = split(",", $header);
shift(@rsids); # Remove 'columns'

#output files
open(PED, ">converted.ped") || die("Could not write to ped file: $!\n");
open(MAP, ">converted.map") || die("Could not write to map file: $!\n");

# write map
foreach my $rsid (@rsids)
{
   if (defined($mapping{$rsid}))
   {
      print MAP join("\t", $mapping{$rsid}{chr}, $rsid, "0", $mapping{$rsid}{pos}) . "\n";
   }
   else
   {
      print MAP join("\t", "0", $rsid, "0", "0") . "\n";
   }
}

# write ped
while (my $line_in = <GENOTYPES>)
{
   chomp $line_in;
   my ($sample_id, @calls) = split(",", $line_in);

   if (defined($metadata{$sample_id}{sex}))
   {
      print PED join("\t", $sample_id, $sample_id, "0", "0", $metadata{$sample_id}{sex}, $metadata{$sample_id}{pheno});
   }
   else
   {
      print PED join("\t", $sample_id, $sample_id, "0", "0", "0", "0");
   }

   my $i = 0;
   foreach my $call (@calls)
   {
      print PED "\t" . conv_call($call, $rsids[$i], \%mapping);
      $i++
   }
   print PED "\n";
}

exit(0);

