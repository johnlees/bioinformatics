#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $missense_regex = qr/missense_variant/;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./vcf_var_counter.pl --vcf <vcf.gz> --lof

Counts damaging variants per sample in a vcf file, annotated with vep

   Options
   --vcf               VCF file to operate on
   --gff               GFF annotation that has been mapped to

   -h, --help          Shows this help.

USAGE

sub print_fasta($$)
{
   my ($id, $trans) = @_;

   open(FASTA, ">$id.fa") || die("Could not write to $id.fa: $!\n");
   print FASTA ">$id\n";
   print FASTA "$trans\n";
   close FASTA;
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($vcf_file, $gff_file, $help);
GetOptions ("vcf=s"  => \$vcf_file,
            "gff=s"  => \$gff_file,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($vcf_file))
{
   print $usage_message;
}
else
{
   # First parse and cache the gff
   open(GFF, $gff_file) || die("Could not open gff file $gff_file: $!\n");
   my %translations;
   while (my $line_in = <GFF>)
   {
      chomp $line_in;

      if ($line_in eq "##FASTA")
      {
         last;
      }
      else
      {
         my @gff_fields = split("\t", $line_in);
         if (defined($gff_fields[2]) && $gff_fields[2] eq "CDS")
         {
            my @info_fields = split(";", $gff_fields[8]);
            my ($id, $translation);
            foreach my $field (@info_fields)
            {
               if ($field =~ /^ID="(.+)"$/)
               {
                  $id = $1;
               }
               elsif ($field =~ /^translation="(.+)"$/)
               {
                  $translation = $1;
               }
            }
            if (defined($id) && defined($translation))
            {
               $translations{$id} = $translation;
            }
         }
      }
   }
   close GFF;

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

         for (my $i = 0; $i < scalar(@csq_fields); $i++)
         {
            if ($csq_fields[$i] =~ /$missense_regex/)
            {
               my @csq_info = split(/\|/, $csq_fields[$i]);

               $csq_info[1] =~ m/^(.+\..+)\..+$/;
               my $id = $1;

               push(@{$vars{$id}{pos}}, $csq_info[7]);
               push(@{$vars{$id}{var}}, $csq_info[8]);
            }
         }
      }
   }
   close VCF;

   # Print info for provean
   foreach my $gene (sort keys %vars)
   {
      if (defined($translations{$gene}))
      {
         print_fasta($gene, $translations{$gene});
      }
      else
      {
         print STDERR "Could not find translation for gene $gene\n";
      }

      open(VARS, ">$gene.vars") || die("Could not write to $gene.vars: $!\n");
      for (my $i = 0; $i < scalar(@{$vars{$gene}{pos}}); $i++)
      {
         if ($vars{$gene}{var}[$i] =~ m/^(.)\/(.)$/)
         {
            my $csq_string = $1 . $vars{$gene}{pos}[$i] . $2;
            print VARS "$csq_string\n";
         }
         else
         {
            print STDERR "Could not parse consequence " . $vars{$gene}{var}[$i] . " at position " . $vars{$gene}{pos}[$i] . "\n";
         }
      }
      close VARS;
   }

}

exit(0);

