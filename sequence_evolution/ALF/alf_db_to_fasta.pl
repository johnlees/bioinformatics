#!/usr/bin/perl -w

# Converts ALF DB sequences into a fasta file for the organism
# Shuffles genes into order (by locus) and removes > lines
# Usage ./alf_db_to_fasta.pl SE0xx_dna.fa > SE0xx.fa

use strict;
use warnings;

use Bio::SeqIO;

my $input_file = $ARGV[0];

open(DB, $input_file) || die("Could not open input $input_file\n");

my @whole_seq;
while (my $db_line = <DB>)
{
   chomp $db_line;
   $db_line =~ m/^>.+, sequence type: type\d+, locus: (\d+)$/;

   my $sequence = <DB>;
   chomp $sequence;
   $whole_seq[$1-1] = $sequence;

   my $blank = <DB>;
}

my $seq_out = Bio::SeqIO->new(-fh=>\*STDOUT,
                              -format=>'Fasta');

$seq_out->write_seq(Bio::Seq->new(-display_id=>($input_file),
                                  -seq=>join("",@whole_seq)));

close DB;

exit(0);

