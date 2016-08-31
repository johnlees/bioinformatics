#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $lof_percentage = 0.8;
my $damaging_regex = qr/stop|frameshift/;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./vcf_var_counter.pl --vcf <vcf.gz> --lof

Counts damaging variants per sample in a vcf file, annotated with vep

   Options
   --vcf               VCF file to operate on

   --lof               If used, require variant to be in the first 80%
                       of the protein
                       CURRENTLY NOT IMPLEMENTED as VEP doesn't annotate with
                       protein length

   -h, --help          Shows this help.

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($vcf_file, $lof, $help);
GetOptions ("vcf=s"  => \$vcf_file,
            "lof"  => \$lof,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($vcf_file))
{
   print $usage_message;
}
else
{
   open(VCF, $vcf_file) || die("Could not open vcf file $vcf_file: $!\n");

   my %vars;
   my @samples;
   while (my $line_in = <VCF>)
   {
      chomp $line_in;
      if ($line_in =~ /^##/)
      {
         next;
      }
      # Read in sample names
      elsif ($line_in =~ /^#/)
      {
         @samples = split("\t", $line_in);
         @samples = splice(@samples, 9);
      }
      else
      {
         # Extract consequences
         my @fields = split("\t", $line_in);
         my @info_fields = split(";", $fields[7]);
         my @csq_fields = split(",", $info_fields[scalar(@info_fields) - 1]);

         my @damaging_fields;
         for (my $i = 0; $i < scalar(@csq_fields); $i++)
         {
            if ($csq_fields[$i] =~ /$damaging_regex/)
            {
               push(@damaging_fields, $i + 1);
            }
         }

         # Count variants with consequences
         for (my $i = 9; $i < scalar(@fields); $i++)
         {
            my @call = split(':', $fields[$i]);
            for (my $j = 0; $j < scalar(@damaging_fields); $j++)
            {
               if ($call[0] eq $damaging_fields[$j])
               {
                  $vars{$samples[$i-9]}++;
                  last;
               }
            }
         }
      }
   }
   foreach my $sample (sort keys %vars)
   {
      print join("\t", $sample, $vars{$sample}) . "\n";
   }

   close VCF;
}

exit(0);

