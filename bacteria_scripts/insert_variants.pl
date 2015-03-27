#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );
use lib dirname( abs_path $0 ) . "/../assembly_scripts";

# Perl modules - assembly
use assembly_common;
use compare_variants;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./insert_variants.pl --map_ref <assembly.fasta> --new_ref <assembly.fasta> --vcf <variants.vcf>

Changes coordinates of vcf between references

   Options
   --map_ref           Reference vcf is mapped to
   --new_ref           Reference to lift over to
   --vcf               VCF file to operate on

   --dirty             Don't clean up temporary files

   -h, --help          Shows this help.

USAGE

# Globals
my $tmp_ref = "reference_renamed.fa";
my $blast_prefix = "blast_windows";

my $new_chrom = "AE007317";

#****************************************************************************************#
#* Functions                                                                            *#
#****************************************************************************************#

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($map_ref, $new_ref, $vcf_in, $dirty, $help);
GetOptions ("map_ref=s"  => \$map_ref,
            "new_ref=s"  => \$new_ref,
            "vcf=s"      => \$vcf_in,
            "dirty"      => \$dirty,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (defined($help))
{
   print STDERR $usage_message;
}
elsif (!defined($map_ref) || !defined($new_ref) || !defined($vcf_in))
{
   print STDERR $usage_message;
}
else
{
   assembly_common::standardise_contig_names($map_ref, $tmp_ref);

   my $blast_output = "$blast_prefix.fa";
   my $variant_lists = compare_variants::extract_vcf_variants($vcf_in);
   compare_variants::variant_windows(100, $variant_lists, $tmp_ref, $blast_output, 1);

   my $blast_scores = compare_variants::blastn_ref($blast_output, $new_ref);

   #copy header then print:
   #vcf lines, tab separated
   #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample
   #AE007317 start-1 . A C . PASS . GT 1
   #
   #so ends up as its own vcf
   #
   system ("bcftools view -h $vcf_in");
   my $num_samples = `bcftools query -f '[%SAMPLE\t]\n' $vcf_in | head -1 | wc -w`;
   chomp $num_samples;

   foreach my $q_id (sort keys %$blast_scores)
   {
      my ($chrom, $pos, $ref, $alt) = split(",", $q_id);

      if ($$blast_scores{$q_id}{start} > $$blast_scores{$q_id}{end})
      {
         $ref = compare_variants::flip_strand($ref, "reverse");
         $alt = compare_variants::flip_strand($alt, "reverse");

         $pos = $$blast_scores{$q_id}{end} - (abs(length($alt) - length($ref)) + 1) - 1;
      }
      else
      {
         $pos = $$blast_scores{$q_id}{start} - 1;
      }

      print join("\t", $new_chrom, $pos, ".", $ref, $alt, ".", "PASS", ".", "GT", "1\t"x$num_samples . "\n");
   }

   unless($dirty)
   {
      unlink $tmp_ref, $blast_output;
   }
}

exit(0);

