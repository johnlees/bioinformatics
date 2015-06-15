#!/usr/bin/perl -w

# Converts ALF DB sequences into a fasta file for the organism
# Shuffles genes into order (by locus) and removes > lines
# Also includes intergenic regions evolved by dawg, and converted with dawg to
# alf
# Usage ./alf_db_to_fasta.pl SE0xx_dna.fa dawg/MSA > SE0xx.fa

use strict;
use warnings;

use Bio::SeqIO;

my $input_file = $ARGV[0];
my $intergenic_folder = $ARGV[1];

open(DB, $input_file) || die("Could not open input $input_file\n");

my @whole_seq;
while (my $db_line = <DB>)
{
   chomp $db_line;
   $db_line =~ m/^>G(\d+)_SE(\d+), sequence type: type\d+, locus: (\d+)$/;
   my $gene_nr = $1;
   my $org_nr = $2;
   my $gene_pos = $3;

   my $sequence = <DB>;
   chomp $sequence;

   if (-e "$intergenic_folder/MSA_i$gene_nr\_dna.fa")
   {
      my $int_in = Bio::SeqIO->new(-file=>"$intergenic_folder/MSA_i$gene_nr\_dna.fa", -format=>'Fasta');
      while (my $int_org_seq = $int_in->next_seq())
      {
         if ($int_org_seq->primary_id() eq "SE$org_nr")
         {
            my $int_region = $int_org_seq->seq();
            $int_region =~ s/-//g;

            $whole_seq[$gene_pos-1] = $int_region;
            last;
         }
      }
   }

   $whole_seq[$gene_pos-1] .= $sequence;

   my $blank = <DB>;
}

my $seq_out = Bio::SeqIO->new(-fh=>\*STDOUT,
                              -format=>'Fasta');

$seq_out->write_seq(Bio::Seq->new(-display_id=>($input_file),
                                  -seq=>join("",@whole_seq)));

close DB;

exit(0);

