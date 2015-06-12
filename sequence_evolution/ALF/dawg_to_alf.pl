#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;

# Converts DAWG output (fasta) into output usable by ALF scripts (MSA and DB)
# usage
# ./dawg_to_alf.pl speciesMapping.txt intergenic_coordinates.txt dawg_output.mfa

my $species_mapping = $ARGV[0];
my $intergenic_coordinates = $ARGV[1];
my $intergenic_results = $ARGV[2];

# Read in name mappings (ALF changes them, DAWG does not)
open(MAP, $species_mapping) || die("Could not open $species_mapping\n");

my %alf_names;
while (my $spec_map_line = <MAP>)
{
   chomp $spec_map_line;

   my ($lane, $alf_name) = split("\t", $spec_map_line);
   $alf_names{$lane} = $alf_name;
}

close MAP;

# Read in coordinates of intergenic regions, and gene they lie before
open(COORDS, $intergenic_coordinates) || die("Could not open $intergenic_coordinates\n");

my %coords;
while (my $coords_line = <COORDS>)
{
   chomp $coords_line;

   my ($gene_nr, $genome_start, $genome_end, $alignment_start, $alignment_end) =
      split("\t", $coords_line);

   $coords{$gene_nr}{start} = $alignment_start;
   $coords{$gene_nr}{end} = $alignment_end;

   my $fasta_out = Bio::SeqIO->new(-file=>">MSA_i$gene_nr\_dna.fa",
                                   -format=>'Fasta');
   $coords{$gene_nr}{file} = $fasta_out;
}

close COORDS;

# Ordered fasta by sample. Take slice for each intergenic region and cat
# together to get mfas for each fasta region, renamed the same as the ALF
# output
my $fasta_in = Bio::SeqIO->new(-file=>$intergenic_results,
                               -format=>'Fasta');
while (my $evolved_seq = $fasta_in->next_seq())
{
   foreach my $gene_nr (keys %coords)
   {
      # great line of code coming up
      $coords{$gene_nr}{file}->write_seq(Bio::Seq(-seq=>$evolved_seq->subseq($coords{$gene_nr}{start}, $coords{$gene_nr}{end}), -display_id=>$alf_names{$evolved_seq->display_id()}));
   }
}

exit(0);

