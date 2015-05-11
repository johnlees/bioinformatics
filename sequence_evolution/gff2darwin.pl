#!/usr/bin/perl -w

#
# Converts gff files to darwin format, for use with ALF
# Uses functions to parse gff v3 from ivr-typer.pl
#

use strict;
use warnings;

# BioPerl modules
use Bio::SeqIO;
use Bio::Tools::GFF;

my $help_message = <<END;
Usage: gff2darwin.pl input.gff > output.db
END

# Need to remove contigs with no features due to failings of
# Bio::Tools::GFF used with attach_seqs
# Should all be sorted by contig:position
sub pre_process_gff($$)
{
   my ($file_in, $file_out) = @_;

   my $fasta_tmp = "tmp.fa";

   if (-e $fasta_tmp)
   {
      die("$fasta_tmp already exists. Will not overwrite\n");
   }

   my $gff_in = Bio::Tools::GFF->new(-file => $file_in,
                                     -gff_version => 3) || die ("Could not open $file_in as gff v3: $!");

   my @contigs;
   my $last_contig = "";

   # Read in all contig sequences, but only add to @contigs if it has at least
   # one feature
   while (my $feature = $gff_in->next_feature())
   {
      if ($feature->seq_id() ne $last_contig)
      {
         $last_contig = $feature->seq_id();
         push(@contigs, $last_contig);
      }
   }

   my @contig_sequences = $gff_in->get_seqs();

   $gff_in->close();

   # Write out contigs with features to a new fasta
   my $fasta_out = Bio::SeqIO->new( -file   => ">$fasta_tmp",
                                    -format => "fasta") || die("Could not write to $fasta_tmp: $!");

   foreach my $contig (reverse @contig_sequences)
   {
      foreach my $included_contig (@contigs)
      {
         if ($contig->display_id() eq $included_contig)
         {
            $fasta_out->write_seq($contig);
         }
      }
   }

   $fasta_out->close();

   # Replace the sequence in the gff. Features and headers remain the same
   #
   # Equivalent command line
   #system("echo \"##FASTA\" | cat $gff_tmp - $fasta_tmp > $file_out");
   open(GFF, "$file_in") || die("Could not open $file_in\n");
   open(FASTA, "$fasta_tmp") || die("Could not open $fasta_tmp\n");
   open(OUT, ">$file_out") || die("Could not write to $file_out\n");

   while (my $gff_line = <GFF>)
   {
      if ($gff_line eq "##FASTA\n")
      {
         last;
      }
      else
      {
         print OUT $gff_line;
      }
   }
   print OUT "\n##FASTA\n";
   while (my $fasta_line = <FASTA>)
   {
      print OUT $fasta_line;
   }

   unlink $fasta_tmp;
   close OUT;

}


#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* Gets input parameters
my $gff_file = $ARGV[0];

if (!defined($gff_file) || !-e $gff_file)
{
   print STDERR "Couldn't find input file\n";
   print STDERR $help_message;
}
else
{
   my $tmp_annot_file = "annotation.tmp";
   pre_process_gff($gff_file, $tmp_annot_file);

   my $gff_in = Bio::Tools::GFF->new(-file => $tmp_annot_file,
                                     -gff_version => 3) || die ("Could not open $tmp_annot_file as gff v3: $!");
   $gff_in->features_attached_to_seqs(1);

   # Look through serially, find any hsd relevant genes
   my @genes;
   while (my $feature = $gff_in->next_feature())
   {
      # Looking at products is more reliable than gene tags
      if ($feature->has_tag("translation"))
      {
         push(@genes, $feature);
      }
   }

   my @sequences = $gff_in->get_seqs();

   $gff_in->close();
   unlink $tmp_annot_file;

   foreach my $gene_feature (@genes)
   {
      my @product = $gene_feature->get_tag_values("product");
      my @trans = $gene_feature->get_tag_values("translation");

      $product[0] =~ s/"//g;
      $trans[0] =~ s/"//g;
      my $sequence = uc($gene_feature->seq->seq());

      my $output_line = "<E><ID>" . $gene_feature->primary_id() . "</ID><DE>" . $product[0]
                        . "</DE><OS>Streptococcus pneumoniae</OS><SEQ>" . $trans[0]
                        . "</SEQ><DNA>" . $sequence . "</DNA></E>";
      print $output_line ."\n";

   }
}

exit(0);
