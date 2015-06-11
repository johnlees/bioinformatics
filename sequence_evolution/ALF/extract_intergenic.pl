#!/usr/bin/perl -w

#
# Extracts intergenic regions from GFF files, for separate evolution by ALF
#

use strict;
use warnings;

# BioPerl modules
use Bio::SeqIO;

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* Gets input parameters
my $gff_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my $output_file = "integenic_regions.fa";

if (!defined($gff_file) || !-e $gff_file || !defined($fasta_file) || !-e $fasta_file)
{
   print STDERR "Couldn't find input files\n";
}
else
{
   my $fasta_in = Bio::SeqIO->new(-file=>$fasta_file,
                              -format=>'Fasta');
   my $fasta_out = Bio::SeqIO->new(-file=>">$output_file",
                              -format=>'Fasta');

   my $in_seq = $fasta_in->next_seq();

   my $region_list = `awk '\$3=="CDS" {print \$0}' $gff_file | grep "translation=" | cut -f 4,5`;
   chomp $region_list;

   my @regions = split("\n", $region_list);
   pop(@regions);

   my $last_end = 0;
   my $concat_pos = 1;
   my $intergenic_seq = "";
   foreach my $gene_region (@regions)
   {
      my ($start, $end) = split("\t", $gene_region);

      if ($start > $last_end+1)
      {
         my $new_intergenic_seq = $in_seq->subseq($last_end+1, $start-1);
         $intergenic_seq .= $new_intergenic_seq;

         my $new_concat_pos = $concat_pos + length($new_intergenic_seq);
         print join("\t", $last_end+1, $start-1, $concat_pos, $new_concat_pos-1) . "\n";
         $concat_pos = $new_concat_pos;
      }

      if ($end > $last_end)
      {
         $last_end = $end;
      }
   }

   $fasta_out->write_seq(Bio::Seq->new(-display_id => 'intergenic_regions', -seq => "$intergenic_seq"));
}

exit(0);
