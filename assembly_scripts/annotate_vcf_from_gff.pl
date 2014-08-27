#!/usr/bin/perl -w
#
use strict;
use warnings;

use Bio::Tools::GFF;

#Globals
my $tmp_annotation = "annotation.tmp";

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
   # Open gff with bioperl, ignoring sequence at the bottom
   my $gff_io = Bio::Tools::GFF->new(-file => $annotations, -gff_version => 3) || die ("Could not open $annotations as gff v3: $!");
   $gff_io->ignore_sequence(1);

   open(ANNOT, ">$tmp_annotation") || die("Could not write to $tmp_annotation: $!\n");
   print ANNOT "CHROM\tFROM\tTO\tTYPE\tID\tGENE\n";
   # Output has tab delimited header (and example lines)
   # CHROM FROM TO TYPE ID GENE
   # contig01 1319 1889 tRNA 2070227_00013 -
   # contig01 1718 2653 CDS 2070227_00018 birA

   while (my $feature = $gff_io->next_feature())
   {
      # First parse contig number
      my $contig;
      if ($feature->seq_id() =~ /^.+\|.+\|(.+)$/)
      {
         $contig = $1;
      }

      #Check if gene is labelled
      my $gene;
      if ($feature->has_tag("gene"))
      {
         my @genes = $feature->get_tag_values("gene");
         $gene = pop(@genes);
      }
      else
      {
         $gene = "";
      }

      # Then print each line
      my $annotate_tmp_line = join("\t", $contig , $feature->start, $feature->end, $feature->primary_tag, $feature->get_tag_values("ID"), $gene);
      print ANNOT "$annotate_tmp_line\n";
   }

   $gff_io->close();
   close ANNOT;

   # Use bcftools to apply this annotation
}

exit(0);

