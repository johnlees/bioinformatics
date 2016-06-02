#!/usr/bin/perl -w

#
# This runs insert_variants on all the intergenic variants (extracted
# previously) and lifts coordinates over to all the S. pneumo carriage isolates
#

use strict;
use warnings;

my @carriage_files = ("/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14355_8_95.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_14.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_20.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_2.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_32.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_38.concat_contig.spades.fa","/lustre/scratch108/bacteria/jl11/nl_carriage/assemblies/14412_3_40.concat_contig.spades.fa");

open(SAMPLES, "intergenic_strep_samples.txt") || die("Could not open intergenic_strep_samples.txt\n");
my @samples;
while (my $line_in = <SAMPLES>)
{
   chomp $line_in;

   push(@samples, $line_in);
}
close SAMPLES;

open(LANES, "strep_pairs.txt") || die("Could not open strep_pairs.txt\n");
my %lanes;
my $header = <LANES>;
while (my $line_in = <LANES>)
{
   chomp $line_in;

   my ($sample, $csf_lane, $blood_lane) = split("\t", $line_in);

   $csf_lane =~ m/^(\d+_\d+)#(\d+)$/;
   $lanes{$sample} = "$1_$2";
}
close LANES;

open(CARRIAGE, "/lustre/scratch108/bacteria/jl11/nl_carriage/carriage_metadata.txt") || die("Could not open /lustre/scratch108/bacteria/jl11/nl_carriage/carriage_metadata.txt\n");
my %carriage_names;
$header = <CARRIAGE>;
while (my $line_in = <CARRIAGE>)
{
   chomp $line_in;

   my ($sample, $sp, $where, $lane) = split("\t", $line_in);

   $sample =~ m/^(\d+)_III$/;;
   $sample = $1;

   $lane =~ m/^(\d+_\d+)#(\d+)$/;
   $lane = "$1_$2";

   $carriage_names{$lane} = $sample;
}
close CARRIAGE;


foreach my $sample (@samples)
{
   foreach my $carriage_assembly (@carriage_files)
   {
      $carriage_assembly =~ m/\/(\d+_\d+_\d+)\.concat_contig\.spades\.fa/;
      my $command = "~/bioinformatics/bacteria_scripts/insert_variants.pl --map_ref /lustre/scratch108/bacteria/jl11/assemblies/$lanes{$sample}/improved_assembly.fa --new_ref $carriage_assembly --vcf $sample/$sample.intergenic.vcf.gz > $sample/$sample.carriage_$carriage_names{$1}.vcf";

      print STDERR $command . "\n";
      system($command);
   }
}

exit(0);

