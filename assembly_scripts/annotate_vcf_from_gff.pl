#!/usr/bin/perl -w
#
use strict;
use warnings;

use vcf_to_gff;

# Headers for bcftools annotation
my $annotation_headers = <<HEADER;
##INFO=<ID=REGION_TYPE,Number=1,Type=String,Description="Type of region variant appears in">
##INFO=<ID=ANNOT_ID,Number=1,Type=String,Description="ID of annotation variant appears in">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name variant appears in">
HEADER

# Usage and input files
my $usage = "./annotate_vcf_from_gff.pl <gff_in> <vcf_file>";

my $annotations = $ARGV[0];
my $vcf_file= $ARGV[1];

# Check input is correct
if (!defined($annotations) || !defined($vcf_file))
{
   print "Usage: $usage\n";
}
elsif (!-e $annotations)
{
   print "Cannot find $annotations\n";
   print "Usage: $usage\n";
}
elsif (!-e $vcf_file)
{
   print "Cannot find $vcf_file\n";
   print "Usage: $usage\n";
}
else
{
   vcf_to_gff::transfer_annotation($annotations, $vcf_file);
}

exit(0);

