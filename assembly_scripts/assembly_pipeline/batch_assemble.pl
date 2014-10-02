#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $wrapper_locations = "/nfs/users/nfs_j/jl11/bioinformatics/assembly_scripts/assembly_pipeline";
my $log_file = "batch_assemble.log";

my $help_message = <<HELP;
Usage: ./batch_assemble.pl <options> --lanes lane_file.txt

Submit assembly and annotation jobs across lsf
Takes a file which contains ids in the form run_lane#tag, one per line

   Options
   --lanes       A list, one per line, of sequencing ids of the form
                 run_lane#tag
   --genus       The genus of the bacteria (used for annotation)

   --output      Output directory (default ./)
   --tmp         Directory to write temporary files to (default /tmp)
   --threads     Number of threads to use (default 1)
   --mem         Comma separated list of memory to submit each step with
                 Default: 4000,2000,1200

   --nosymlinks  Do not create subdirectories with symlinks, as they
                 already exist
   --improve     Run from improvement step
   --annotate    Run from annotation step

   --help      Displays this help message

HELP

# bsub parameters
my $spades_threads = 1;
my $spades_mem = 4000;
my $spades_tmp = 2000;

my $improve_mem = 2000;
my $improve_tmp = 1300;

my $annotate_threads = 1;
my $annotate_mem = 1200;

#**********************************************************************#
#* Subs                                                               *#
#**********************************************************************#
sub get_fastq($$)
{
   my ($lane, $output_dir) = @_;

   my $fastq_dir = `pathfind -t lane -i $lane`;
   chomp($fastq_dir);

   symlink "$fastq_dir/$lane" . "_1.fastq.gz", "$output_dir/$output_dir" . "_1.fastq.gz";
   symlink "$fastq_dir/$lane" . "_2.fastq.gz", "$output_dir/$output_dir" . "_2.fastq.gz";
}

sub run_getid($)
{
   my ($command) = @_;

   my $submit = `$command`;

   my $jobid;
   if ($submit =~ /^Job <(\d+)>/)
   {
      $jobid = $1;
   }
   else
   {
      print STDERR "Command $command could not be submitted\n";
      $jobid = 0;
   }

   return($jobid);
}

#**********************************************************************#
#* Main                                                               *#
#**********************************************************************#

# Get input parameters
my ($sample_file, $genus, $output_directory, $threads, $mem_string, $from_improvement, $from_annotation, $no_symlinks, $tmp_directory, $help);
GetOptions("lanes=s" => \$sample_file,
           "output|o=s" => \$output_directory,
           "genus=s" => \$genus,
           "mem=s" => \$mem_string,
           "threads=i" =>\$threads,
           "tmp=s" => \$tmp_directory,
           "nosymlinks" => \$no_symlinks,
           "improve" => \$from_improvement,
           "annotate" => \$from_annotation,
           "help" => \$help) || die("$!\n$help_message");

if (!defined($sample_file) || !-e $sample_file)
{
   print STDERR "Input sample file not given, or does not exist\n";
   print STDERR $help_message;
}
elsif (defined($help))
{
   print $help_message;
}
else
{
   # Parse input
   if (!defined($output_directory))
   {
      $output_directory = "./";
   }

   if (!defined($tmp_directory))
   {
      $tmp_directory = "/tmp/spades";
   }

   if (defined($threads))
   {
      $spades_threads = $threads;
      $annotate_threads = $threads;
   }

   if (defined($mem_string))
   {
      ($spades_mem, $improve_mem, $annotate_mem) = split(",", $mem_string);
   }

   # Set up input files
   my @samples;
   open(LANES, "$sample_file") || die("Could not open $sample_file");

   while (my $lane = <LANES>)
   {
      chomp($lane);

      if ($lane =~ m/^(\d+)_(\d+)[#_](\d+)$/)
      {
         my $lane_name = "$1_$2_$3";
         push(@samples, $lane_name);

         # Rename to avoid irritating hashes, then symlink fastqs into that
         # directory
         if (!-d $lane_name)
         {
            mkdir $lane_name;
            mkdir "$lane_name/logs";
         }

         unless($no_symlinks)
         {
            get_fastq($lane, $lane_name);
         }
      }
      else
      {
         print STDERR "$lane doesn't look like run_lane#tag\nSkipping\n";
      }
   }

   close LANES;

   open(LOG, ">$log_file") || die("Couldn't write to log file $log_file: $!\n");
   print LOG "Lane\tAssembly_job\tImprovement_job\tAnnotation_job\n";

   my $spades_jobid = "na";
   my $improve_jobid = "na";
   my $annotation_jobid = "na";

   foreach my $sample (@samples)
   {
      my $forward_reads = "$sample/$sample" . "_1.fastq.gz";
      my $reverse_reads = "$sample/$sample" . "_2.fastq.gz";

      my $tmp_sample_dir = "$tmp_directory/$sample";

      # Submit assembly job & filter
      my $assembly_command;
      if ($tmp_directory =~ /^\/tmp/)
      {
         $assembly_command = "bsub -o $sample/logs/spades.%J.o -e $sample/logs/spades.%J.e -n$spades_threads -R \"span[hosts=1]\" -R \"select[mem>$spades_mem] rusage[mem=$spades_mem]\" -R \"select[tmp>$spades_tmp]\" -M$spades_mem $wrapper_locations/spades_wrapper.pl $forward_reads $reverse_reads $sample $spades_threads $tmp_sample_dir";
      }
      else
      {
         $assembly_command = "bsub -o $sample/logs/spades.%J.o -e $sample/logs/spades.%J.e -n$spades_threads -R \"span[hosts=1]\" -R \"select[mem>$spades_mem] rusage[mem=$spades_mem]\" -M$spades_mem $wrapper_locations/spades_wrapper.pl $forward_reads $reverse_reads $sample $spades_threads $tmp_sample_dir";
      }

      if(!defined($from_annotation) && !defined($from_annotation))
      {
         $spades_jobid = run_getid($assembly_command);
      }

      # Improvement contingent on success
      my $improvement_bsub = "bsub -o $sample/logs/improve_assembly.%J.o -e $sample/logs/improve_assembly.%J.e -R \"select[mem>$improve_mem] rusage[mem=$improve_mem]\" -M$improve_mem";
      if ($spades_jobid ne "na")
      {
         $improvement_bsub .= " -w \"done($spades_jobid)\"";
      }
      if ($tmp_directory =~ /^\/tmp/)
      {
         $improvement_bsub .= " -R \"select[tmp>$improve_tmp]\"";
      }

      my $improvement_command = "$wrapper_locations/improvement_wrapper.pl $forward_reads $reverse_reads $sample/scaffolds.filtered.fasta $sample";

      if (!defined($from_annotation))
      {
         $improve_jobid = run_getid("$improvement_bsub $improvement_command");
      }

      # Annotation contingent on success
      my $annotation_bsub = "cd $sample && bsub -o logs/annotate.%J.o -e logs/annotate.%J.e -R \"select[mem>$annotate_mem] rusage[mem=$annotate_mem]\" -M$annotate_mem -n$annotate_threads -R \"span[hosts=1]\"";
      if ($improve_jobid ne "na")
      {
         $annotation_bsub .= "-w \"done($improve_jobid)\"";
      }

      my $annotation_command = "$wrapper_locations/annotation_wrapper.pl improved_assembly.fa $sample $genus $annotate_threads";
      my $annotation_jobid = run_getid("$annotation_bsub $annotation_command");

      my $report = join("\t", $sample, $spades_jobid, $improve_jobid, $annotation_jobid);

      print LOG "$report\n";
   }

   close LOG;
   print STDERR "Job details are in $log_file\n";
}

exit(0);

