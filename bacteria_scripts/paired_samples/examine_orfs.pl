#!/usr/bin/perl -w

# Takes ORFs from summarise_mutation.pl, and creates a multi-fasta of them
# for analysis by blast

use strict;
use warnings;

use Bio::SeqIO;

my $annotation_root = "/lustre/scratch108/bacteria/jl11/assemblies";

my $orf_file = $ARGV[0];
open(ORFS, $orf_file) || die("Couldn't open $orf_file\n");

my $sequence_out = Bio::SeqIO->new( -file   => ">orfs.fa",
                                    -format => "fasta") || die ($!);

while (my $orf_line = <ORFS>)
{
   chomp $orf_line;

   my ($gene_name, $count) = split("\t", $orf_line);
   if ($gene_name =~ /^ORF_(\d+_\d+_\d+)_(\d+)/)
   {
      my $lane = $1;
      my $orf = "$1_$2";

      my $annotation = "$annotation_root/$lane/annotation/$lane.ffn";

      my $sequence_in = Bio::SeqIO->new( -file   => "<$annotation",
                                      -format => "fasta" ) || die ($!);

      while (my $sequence = $sequence_in->next_seq())
      {
         if ($orf eq $sequence->display_id())
         {
            $sequence_out->write_seq($sequence);
            last;
         }
      }
   }
}

close ORFS;

exit(0);
