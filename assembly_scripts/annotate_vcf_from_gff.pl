#!/usr/bin/perl -w
#
use strict;
use warnings;

use Bio::Tools::GFF;

#Globals
my $tmp_annotation = "annotation.tmp";
my $annotation_header_file = "annotation_headers.tmp";

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
   # Open gff with bioperl, ignoring sequence at the bottom
   # EDIT: Ignoring sequence seems to cause loop to hang, so include it for now
   # despite overheads
   my $gff_io = Bio::Tools::GFF->new(-file => $annotations, -gff_version => 3) || die ("Could not open $annotations as gff v3: $!");
   $gff_io->ignore_sequence(0);

   open(ANNOT, ">$tmp_annotation") || die("Could not write to $tmp_annotation: $!\n");
   print ANNOT "CHROM\tFROM\tTO\tREGION_TYPE\tANNOT_ID\tGENE\n";
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
   # First add the necessary new headers
   open (HEADERS, ">$annotation_header_file") || die("Could not write to $annotation_header_file: $!\n");
   print HEADERS $annotation_headers;
   close HEADERS;

   # bgzip and tabix annotations file
   system("bgzip $tmp_annotation");
   system("tabix -s 1 -b 2 -e 3 $tmp_annotation.gz");

   # Also annotate allele counts at the same time
   my $bcftools_command = "bcftools annotate -p fill-AN-AC -a $tmp_annotation.gz -h $annotation_header_file -c CHROM,FROM,TO,REGION_TYPE,ANNOT_ID,GENE -o $vcf_file -O z $vcf_file";
   system($bcftools_command);
   system("bcftools index $vcf_file");

   # Finally, remove temporary files
   unlink "$tmp_annotation.gz", "$tmp_annotation.gz.tbi", $annotation_header_file;
}

exit(0);

