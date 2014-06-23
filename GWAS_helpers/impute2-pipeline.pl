#!/usr/bin/perl -w

use strict;
use warnings;

use POSIX;

#
# Globals
#
my $wait_time = 30;
my $default_memory = 9000;
my $mem_increment = 1500;
my $chunk_length = 5000000;

my $job_id_file = "job_ids.log";

# file locations
my $input_prefix = "shapeit2/cases-ALS-BPROOF.shapeit2.chr";
my $ref_directory = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono";
my $output_directory = "impute2";

my $map_prefix = "genetic_map_chr";
my $map_suffix = "_combined_b37.txt";

my $ref_prefix = "ALL.chr";
my $hap_suffix = ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz";
my $leg_suffix = ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz";
my $ref_sample_file = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample";

my $X_prefix = "ALL_1000G_phase1integrated_v3_chrX_nonPAR";
my $X_hap_suffix = "_impute.hap.gz";
my $X_leg_suffix = "_impute.legend.gz";

#****************************************************************************************#
#* Subs                                                                                 *#
#****************************************************************************************#
sub chrom_jobs($)
{
   my ($chrom_file) = @_;

   my $last_line;

   open(LEGEND, "gzip -dc $chrom_file |") or die("Could not open $chrom_file\n");
   while (my $legend_line = <LEGEND>)
   {
      # Get to end of file
      $last_line = $legend_line if eof;
   }

   my ($rsid, $position, $a0, $a1, $var_type, $seq_type, @pop_maf) = split(/\s+/, $last_line);
   my $num_jobs = ceil($position / $chunk_length);

   close LEGEND;

   return($num_jobs);

}

#
# job finished/memory overrun function
#
sub check_status($)
{
   my ($jobid) = @_;

   my $status;

   my $bjobs = `bjobs -a -noheader -o "stat exit_code delimiter=','" $jobid`;
   my ($bjobs_stat, $exit_code) = split(/,/, $bjobs);

   if ($bjobs_stat eq "RUN" || $bjobs_stat eq "PEND")
   {
      $status = "RUN";
   }
   elsif ($bjobs_stat eq "EXIT" && $exit_code == 130)
   {
      my $mem_over = `bjobs -l $jobid | grep -l "TERM_MEMLIMIT" | wc -l`;
      chomp($mem_over);

      if ($mem_over == 1)
      {
         $status = "MEMLIMIT";
      }
      else
      {
         $status = "Status: $bjobs_stat. Exit code: $exit_code. NOT memlimit exceeded";
      }
   }
   elsif ($bjobs_stat eq "DONE" && $exit_code == 0)
   {
      $status = "DONE";
   }
   else
   {
      $status = "Status: $bjobs_stat. Exit code: $exit_code";
   }

   return $status;
}

# impute2 command
sub run_impute2($$$)
{
   my ($chr, $chunk, $memory) = @_;

   my $jobid;

   # e.g. ./impute2 \
   # -m ./Example/example.chr22.map \
   # -h ./Example/example.chr22.1kG.haps \
   # -l ./Example/example.chr22.1kG.legend \
   # -g ./Example/example.chr22.study.gens \
   # -strand_g ./Example/example.chr22.study.strand \
   # -int 20.4e6 20.5e6 \
   # -Ne 20000 \
   # -o ./Example/example.chr22.one.phased.impute2
   my $int_start = ($chunk - 1) * $chunk_length;
   my $int_end = $int_start + $chunk_length;

   my ($m_file, $g_file, $g_sample_file, $h_file, $l_file, $impute2_command);
   unless ($chr eq "X")
   {
      $h_file = "$ref_directory/$ref_prefix$chr$hap_suffix";
      $l_file = "$ref_directory/$ref_prefix$chr$leg_suffix";
      $m_file = "$ref_directory/$map_prefix$chr$map_suffix";
      $g_file = "$input_prefix$chr.haps";
      $g_sample_file = "$input_prefix$chr.sample";

      $impute2_command = "impute2 -m $m_file -h $h_file -l $l_file -known_haps_g $g_file -sample_g $g_sample_file -int $int_start $int_end -Ne 20000 -o $output_directory/$ref_prefix.impute2.$chr.$chunk";
   }
   else
   {
      $h_file = "$ref_directory/$X_prefix$chr$X_hap_suffix";
      $l_file = "$ref_directory/$X_prefix$chr$X_leg_suffix";
      $m_file = "$ref_directory/$map_prefix" . "X_nonPAR" . "$map_suffix";
      $g_file = "$input_prefix" . "X.haps";
      $g_sample_file = "$input_prefix" . "X.sample";

      $impute2_command = "impute2 -chrX -m $m_file -h $h_file -l $l_file -known_haps_g $g_file -sample_g $g_sample_file -int $int_start $int_end -Ne 20000 -o $output_directory/$ref_prefix.impute2.$chr.$chunk";
   }

   my $bsub_command = "bsub -J " . '"impute2"' . " -o $output_directory/impute2.%J.$chr.$chunk.o -e $output_directory/impute2.%J.$chr.$chunk.e -R " . '"' . "select[mem>$memory] rusage[mem=$memory]" . '"' . " -M$memory -q long";

   # example output: Job <5521290> is submitted to queue <normal>.
   my $submit = `$bsub_command $impute2_command`;

   if ($submit =~ /^Job <(\d+)>/)
   {
      $jobid = $1;
   }
   else
   {
      $jobid = 0;
   }

   return($jobid);
}

sub the_time()
{
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
   $year = sprintf("%02d", $year % 100);
   $mon += 1;

   my $time_string = "$year-$mon-$mday/$hour:$min:$sec";

   return $time_string;
}

sub update_job_id_file($)
{
   my ($job_hash_ref) = @_;
   my %job_hash = %$job_hash_ref;

   open(JOBS, ">$job_id_file") || die("Could not open $job_id_file: $!");

   foreach my $chr (sort keys %job_hash)
   {
      foreach my $chunk (sort keys %{$job_hash{$chr}})
      {
         my $job_line = join("\t", $job_hash{$chr}{$chunk}, $chr, $chunk);
         print JOBS "$job_line\n";
      }
   }

   close JOBS;
}

sub read_job_id_file()
{
   my %job_ids;
   open(JOBS, $job_id_file) || die("Could not open $job_id_file: $!");

   while (my $job_line = <JOBS>)
   {
      chomp($job_line);

      my ($job_id, $chr, $chunk) = split(/\t/, $job_line);
      $job_ids{$chr}{$chunk} = $job_id;
   }

   close JOBS;

   return(\%job_ids);
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

open(LOG, ">>impute2-pipeline.log") || die ("Could not open impute2-pipeline.log for writing");
open(ERRORS, ">>impute2-pipeline.err") || die ("Could not open impute2-pipeline.err for writing");

# Get number of jobs per chromosome
my %num_jobs;

for (my $i = 1; $i <= 22; $i++)
{
   my $file_name = "$ref_directory/$ref_prefix$i$leg_suffix";
   $num_jobs{$i} = chrom_jobs($file_name);
   print LOG "Chr $i has $num_jobs{$i} chunks of 5Mb\n";
   print "Chr $i has $num_jobs{$i} chunks of 5Mb\n";
}
my $x_chr_name = "$ref_directory/$X_prefix$X_leg_suffix";
$num_jobs{"X"} = chrom_jobs($x_chr_name);
print LOG "Chr X has " . $num_jobs{"X"} . " chunks of 5Mb\n";
print "Chr X has " . $num_jobs{"X"} . " chunks of 5Mb\n";

my %jobid;
my %memory;

# Set up jobs to run, then loop through them
my @chrs = (1..22);
push(@chrs, "X");

if (!-e $job_id_file)
{
   foreach my $chr (@chrs)
   {
      for (my $chunk = 1; $chunk <= $num_jobs{$chr}; $chunk++)
      {
         $jobid{$chr}{$chunk} = 0;
      }
   }
}
else
{
   my $jobs_in = read_job_id_file();
   %jobid = %$jobs_in;
}

update_job_id_file(\%jobid);

my $imputation_ongoing = 1;
my %processed;

while ($imputation_ongoing)
{
   $imputation_ongoing = 0;

   foreach my $chr (sort keys %jobid)
   {
      foreach my $chunk (sort keys %{$jobid{$chr}})
      {
         if ($jobid{$chr}{$chunk} == 0)
         {
            $memory{$chr}{$chunk} = $default_memory;
            $jobid{$chr}{$chunk} = run_impute2($chr, $chunk, $memory{$chr}{$chunk});
            update_job_id_file(\%jobid);
            $imputation_ongoing = 1;
         }
         else
         {
            if (!defined($processed{$chr}{$chunk}))
            {
               my $status = check_status($jobid{$chr}{$chunk});
               if ($status eq "RUN")
               {
                  $imputation_ongoing = 1;
               }
               elsif ($status eq "MEMLIMIT")
               {
                  if (!defined($memory{$chr}{$chunk}) || $memory{$chr}{$chunk} <= $default_memory)
                  {
                     $memory{$chr}{$chunk} = $mem_increment + $default_memory;
                  }
                  else
                  {
                     $memory{$chr}{$chunk} += $mem_increment;
                  }
                  print LOG "Chromosome:$chr chunk:$chunk ran over memory limit, resubmitting with " . $memory{$chr}{$chunk} . "MB\n";
                  print "Chromosome:$chr chunk:$chunk ran over memory limit, resubmitting with " . $memory{$chr}{$chunk} . "MB\n";
                  $jobid{$chr}{$chunk} = run_impute2($chr, $chunk, $memory{$chr}{$chunk});
                  update_job_id_file(\%jobid);
                  $imputation_ongoing = 1;
               }
               elsif ($status eq "DONE")
               {
                  print LOG "Chromosome:$chr chunk:$chunk complete at " . the_time() . "\n";
                  print "Chromosome:$chr chunk:$chunk complete at " . the_time() . "\n";
                  $processed{$chr}{$chunk} = 1;
               }
               else
               {
                  print ERRORS "Chromosome:$chr chunk:$chunk jobid $jobid{$chr}{$chunk} failed with unknown error: $status\n";
                  print "Chromosome:$chr chunk:$chunk jobid $jobid{$chr}{$chunk} failed with unknown error: $status\n";
                  $processed{$chr}{$chunk} = 1;
               }
            }
         }
      }
   }
   print LOG "Jobs checked. Waiting $wait_time seconds\n\n";
   print "Jobs checked. Waiting $wait_time seconds\n\n";
   sleep($wait_time);
}

print LOG "All jobs finished - creating final output \n\n";
print "All jobs finished - creating final output \n\n";

# Concat output files

foreach my $chr (sort keys %jobid)
{
   my $cat_command = "cat ";
   for (my $chunk = 1; $chunk <= $num_jobs{$chr}; $chunk++)
   {
      $cat_command .= "$output_directory/$ref_prefix.impute2.$chr.$chunk.gen ";
   }
   my $output_file = "$output_directory/$ref_prefix.impute2.chr$chr.gen";
   $cat_command .= "> $output_file";

   my $cat_result = `$cat_command`;
   print LOG "Chromosome $chr output: $output_file result: $cat_result\n";
   system("gzip $output_file");
}

close LOG;
close ERRORS;

exit(0);

