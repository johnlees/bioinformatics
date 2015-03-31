#!/usr/bin/perl -w

use strict;
use warnings;

# USAGE
# ./diff_vep_batch.pl > diff_consequences.txt 2> diff_vep_batch.log

my $root_dir = "/lustre/scratch108/bacteria/jl11/";

open(SAMPLES, "samples.txt") || die("Could not open samples.txt\n");

while (my $sample = <SAMPLES>)
{
   chomp($sample);
   chdir($sample);
   print STDERR "Sample $sample\n";

   # get ref mapped to
   open(READS, "$root_dir/pairs_analysis/strep_pneumo_pairs/$sample/reads.txt") || die("Could not open reads.txt for $sample\n");

   my $ref_line = <READS>;
   chomp($ref_line);
   my ($sample_tissue, $lane_1, $lane_2) = split("\t", $ref_line);

   $lane_1 =~ m/\/(\d+)_(\d+)#(\d+)_1.fastq.gz$/;
   my $map_ref = "$root_dir/assemblies/$1_$2_$3/improved_assembly.fa";

   close READS;

   # Transfer to new ref
   system("perl ~/bioinformatics/bacteria_scripts/insert_variants.pl --vcf ../../strep_pneumo_pairs/$sample/$sample.vcf.gz --new_ref ../Streptococcus_pneumoniae_R6_v1.fa --map_ref $map_ref > diffs.vcf");

   my %vars;
   my $var_output = `bcftools view -H diffs.vcf`;
   my @var_lines = split("\n", $var_output);

   foreach my $var_line (@var_lines)
   {
      my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $csf, $blood) = split("\t", $var_line);

      if ($csf == 1)
      {
         $vars{"$chr:$pos"} = "csf";
      }
      else
      {
         $vars{"$chr:$pos"} = "blood";
      }
   }

   # Run vep
   if (scalar(keys %vars) != 0)
   {
      system("perl /nfs/users/nfs_j/jl11/installations/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl -i diffs.vcf --cache --species Streptococcus_pneumoniae_R6 --cache_version 25 --offline --output_file diff_effects --coding_only --force_overwrite &> vep_diff.log");

      # Summarise output
      open(VEP, "diff_effects") || die("Couldn't open vep output for sample $sample\n");
      while (my $line_in = <VEP>)
      {
         chomp $line_in;

         unless ($line_in =~ /^#/)
         {
            my ($var, $loc, $allele, $gene, $feature, $feature_type, $consequence, $cDNA_position, $CDS_position, $Protein_position, $Amino_acids, $Codons, $Existing_variation, $Extra) = split("\t", $line_in);

            if ($feature =~ /dlt/)
            {
               $loc =~ m/^(.+):(\d+)/;
               print join("\t", $sample, $feature, $consequence, $vars{"$1:$2"}) . "\n";
            }
         }
      }
   }

   chdir("..");
}

close SAMPLES;

exit(0);

