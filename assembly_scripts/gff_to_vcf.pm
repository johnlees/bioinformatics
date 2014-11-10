#!/usr/bin/perl -w
#
package gff_to_vcf;

use strict;
use warnings;

use Bio::Tools::GFF;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use assembly_common;

#Globals
my $tmp_annotation = "annotation.tmp";
my $annotation_header_file = "annotation_headers.tmp";

# Global locations of software needed
my $bcftools = "/nfs/users/nfs_j/jl11/software/bin/bcftools";
my $tabix_location = "/usr/bin/tabix"; # Must be <v1.0!
my $bgzip_location = "bgzip";

# Headers for bcftools annotation
my $annotation_headers = <<HEADER;
##INFO=<ID=REGION_TYPE,Number=1,Type=String,Description="Type of region variant appears in">
##INFO=<ID=ANNOT_ID,Number=1,Type=String,Description="ID of annotation variant appears in">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name variant appears in">
HEADER

sub transfer_annotation($$)
{
   my ($annotations, $vcf_file) = @_;

   my $stderr_file = "bcf_annotate.err";

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
   system("$bgzip_location $tmp_annotation");
   system("$tabix_location -s 1 -b 2 -e 3 $tmp_annotation.gz");

   # Also annotate allele counts at the same time
   # bcftools v1.1: AN, AC already calculated. Plugins available from bcftools
   # plugin
   my $bcftools_command = "$bcftools annotate -a $tmp_annotation.gz -h $annotation_header_file -c CHROM,FROM,TO,REGION_TYPE,ANNOT_ID,GENE -o $vcf_file -O z $vcf_file 2> $stderr_file";
   system($bcftools_command);

   # Need to refresh index
   system("$bcftools index -f $vcf_file 2>> $stderr_file");

   # Finally, remove temporary files
   unlink "$tmp_annotation.gz", "$tmp_annotation.gz.tbi", $annotation_header_file;
   assembly_common::add_tmp_file($stderr_file);
}

# Creates a file of exons for use with frameshift annotation
sub create_exons_tab($$)
{
   my ($gff_file, $exons_file) = @_;

   # Input and output
   open (GFF, "$gff_file") || die("Could not open $gff_file\n");
   open (EXONS, ">$exons_file") || die("Could not open $exons_file for writing\n");

   # Headers in output are CHROM, FROM, TO
   while (my $gff_line = <GFF>)
   {
      chomp $gff_line;

      # Reached end of annotation?
      if ($gff_line eq "##FASTA")
      {
         last;
      }
      elsif ($gff_line =~ /^##/)
      {
         next;
      }
      else
      {
         my ($contig, $predictor, $region, $start, $end, @gff_fields) = split("\t", $gff_line);

         if ($region eq "CDS")
         {
            $contig =~ m/^.+\|SC\|(.+)$/;
            my $chrom = $1;

            print EXONS join("\t", $chrom, $start, $end) . "\n";
         }

      }
   }

   close GFF;
   close EXONS;

   # Output is bgzip compressed and tabix indexed
   system("bgzip $exons_file");
   system("tabix -s1 -b2 -e3 $exons_file.gz");
}

1;

