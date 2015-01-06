#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

#
# Given a list of illumina lanes, runs both the mapping and assembly
# based reference free variant callers and produces a final call set
# for each sample (i.e. two sequences - in my case blood and csf
# isolates)
#

my $usage_message = <<USAGE;
Usage: ./pairs_pipe.pl --lanes <list_of_lanes.txt> --assembly_dir <assemblies>

Uses cortex to call indels and samtools/bcftools to call snps between a pair
of samples

   Options
   (Required)
   --lanes           A list of lanes to run the pipeline on. Three tab
                     separated columns: sample ID, lane 1, lane 2
   --assembly_dir    A directory containing assemblies produced by
                     batch_assemble.pl for all lanes in the file

   (Optional)
   --output          Output directory. Default ./

   -h, --help        Shows this help.

USAGE

my $cortex_mem = 2000;
my $map_memory = 4000;

my $job_regex = qr/^Job <(\d+)>/;

#**********************************************************************#
#* Subs                                                               *#
#**********************************************************************#
sub get_fastq($$)
{
   my ($lane, $output_dir) = @_;

   my $fastq_dir = `pathfind -t lane -i $lane`;
   chomp($fastq_dir);

   symlink "$fastq_dir/$lane" . "_1.fastq.gz", "$output_dir/$lane" . "_1.fastq.gz";
   symlink "$fastq_dir/$lane" . "_2.fastq.gz", "$output_dir/$lane" . "_2.fastq.gz";
}

#**********************************************************************#
#* Main                                                               *#
#**********************************************************************#
my ($lane_file, $assembly_directory, $output_directory, $help);
GetOptions("lanes=s" => \$lane_file,
           "output|o=s" => \$output_directory,
           "assembly_dir=s" => \$assembly_directory,
           "help|h" => \$help) || die("$!\n$usage_message");

# Parse options
if (defined($help) || !defined($lane_file) || !defined($assembly_directory))
{
   print STDERR $usage_message;
}
else
{
   if (!defined($output_directory))
   {
      $output_directory = "./";
   }
   chdir $output_directory;

   open(JOBS, ">pairs_pipe.txt") || die("Could not write to pairs_pipe.txt\n");
   print JOBS join("\t", "Sample", "Cortex ID", "Map ID", "Concat ID") . "\n";

   open(LANES, $lane_file) || die("Could not open $lane_file\n");

   while (my $lane_line = <LANES>)
   {
      chomp $lane_line;
      my ($sample, $csf_lane, $blood_lane) = split("\t", $lane_line);

      # Produce directories for each sample
      if (!-d $sample)
      {
         mkdir $sample;
         get_fastq($csf_lane, $sample);
         get_fastq($blood_lane, $sample);

         chdir $sample;

         open(READS, ">reads.txt") || die("Could not write to $sample/reads.txt\n");
         print READS join("\t", $sample . "_csf", "$csf_lane" . "_1.fastq.gz", "$csf_lane" . "_2.fastq.gz\n" .
                                $sample . "_blood", "$blood_lane" . "_1.fastq.gz", "$blood_lane" . "_2.fastq.gz");
         close READS;
      }
      else
      {
         chdir $sample;
      }

      if (!-d "logs")
      {
         mkdir "logs";
      }

      # bsub cortex
      mkdir "cortex";
      chdir "cortex";

      my $bsub_cortex = "bsub -o ../logs/cortex.%J.o -e ../logs/cortex.%J.e -R \"select[mem>$cortex_mem] rusage[mem=$cortex_mem]\" -M$cortex_mem";
      my $cortex_command = "~/bioinformatics/assembly_scripts/reference_free_variant_caller.pl --cortex -a $assembly_directory/$csf_lane/improved_assembly.fa -g $assembly_directory/$csf_lane/annotation/annotation.gff --separate-correct -r reads.txt -o $sample";

      my $cortex_job = `$bsub_cortex $cortex_command`;
      # Job <3849944> is submitted to default queue <normal>
      $cortex_job =~ $job_regex;
      my $job1_id = $1;

      chdir "..";

      # bsub mapping
      mkdir "mapping";
      chdir "mapping";

      my $bsub_mapping = "bsub -o ../logs/mapping.%J.o -e ../logs/mapping.%J.e -R \"select[mem>$map_memory] rusage[mem=$map_memory]\" -M$map_memory";
      my $map_command = "perl ~/bioinformatics/bacteria_scripts/map_snp_call.pl -a -a $assembly_directory/$csf_lane/improved_assembly.fa -g $assembly_directory/$csf_lane/annotation/annotation.gff -r reads.txt -o $sample -p 1e-6";
      my $mapping_job = `$bsub_mapping $map_command`;

      $mapping_job =~ $job_regex;
      my $job2_id = $1;

      chdir "..";

      # bsub concat command, conditional on jobs finishing
      my $concat_command = "bsub -w \"done($job1_id)\" -w \"done($job2_id)\" -o logs/combine.%J.o ~/bioinformatics/bacteria_scripts/paired_samples/combine_output.pl";
      my $concat_job = `$concat_command`;

      $concat_job =~ $job_regex;
      my $job3_id = $1;

      print JOBS join("\t", $sample, $job1_id, $job2_id, "$job3_id\n");

      chdir "..";
   }
   close LANES;
   close JOBS;

   print STDERR "Job details written to pairs_pipe.txt\n"

}

exit(0);

