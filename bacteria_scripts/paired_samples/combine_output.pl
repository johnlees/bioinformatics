#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = "./combine_output.pl --cortex sample/cortex --map sample/mapping --output sample";

my ($cortex_dir, $map_dir, $output_prefix, $help);
GetOptions("cortex=s" => \$cortex_dir,
           "map=s" => \$map_dir,
           "output=s" => \$output_prefix,
           "help|h" => \$help) || die("$!\n$usage_message");

if (defined($help))
{
   print STDERR $usage_message;
}
else
{
   system("bcftools view -v indels $cortex_dir/$output_prefix.filtered_calls.vcf.gz -O z -o $output_prefix/$output_prefix.indels.vcf.gz");
   system("bcftools view -v snps $map_dir/$output_prefix.diff.vcf.gz -O z -o $output_prefix/$output_prefix.snps.vcf.gz");

   system("bcftools concat - a$output_prefix/$output_prefix.snps.vcf.gz $output_prefix/$output_prefix.indels.vcf.gz -O z -o $output_prefix/$output_prefix.vcf.gz");
   system("bcftools index $output_prefix/$output_prefix.vcf.gz");

   unlink "$output_prefix/$output_prefix.indels.vcf.gz", "$output_prefix/$output_prefix.snps.vcf.gz";
}


exit(0);

