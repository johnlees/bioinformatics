#!/usr/bin/perl -w

use strict;
use warnings;

my $jobs = 300;
my $job_nr = $ENV{'LSB_JOBINDEX'};

for (my $i = ($job_nr-1)*$jobs+1; $i <= $jobs*$job_nr; $i++)
{
   if ($i > 3069)
   {
      last;
   }

   my $sed_command = "sed '$i\q;d' /lustre/scratch108/bacteria/jl11/seer/maela_kmers/seer_new/assemblies.txt";
   my $assembly_file = `$sed_command`;
   chomp $assembly_file;

   $assembly_file =~ m/(\d+_\d+#\d+)\/velvet_assembly\/contigs\.fa$/;
   my $assembly_base = $1;

   system("cp $assembly_file $assembly_base.fa");
   system("bwa index $assembly_base.fa");
   system("bwa fastmap -l 15 $assembly_base.fa trimethoprim_hits.fa > $assembly_base.smem");
   system("./smem_to_kmer.pl trimethoprim_hits.fa $assembly_base.smem > $assembly_base.txt");

   unlink glob "$assembly_base.fa*";
}

exit(0);

