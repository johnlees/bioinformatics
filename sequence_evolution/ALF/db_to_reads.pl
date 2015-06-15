#!/usr/bin/perl -w

# Converts ALF DB files into reads using pIRS, then assembles them into contigs using
# velvet optimiser
#
# Run
# bsub -o logs/assemble.%J.%I.o -e logs/assemble.%J.%I.e -R "select[mem>500] rusage[mem=500]" -M500 -J "reads[1-3069]%1000" ./db_to_reads.pl

use strict;
use warnings;

use File::Path qw(remove_tree);

my $job_id = $ENV{'LSB_JOBINDEX'};
if (!defined($job_id))
{
   $job_id = $ARGV[0];
}

# Convert DB to fasta
my $file_name;
if ($job_id < 1000)
{
   $file_name = sprintf("SE%03d", $job_id);
}
else
{
   $file_name = "SE$job_id\_dna.fa";
}
system("~/bioinformatics/sequence_evolution/ALF/alf_db_to_fasta.pl DB/$file_name\_dna.fa > DB/$file_name\_genome.fa");

# Generate reads
system("~/installations/pIRS_111/pirs simulate -i DB/$file_name\_genome.fa -l 100 -x 50 -m 250 -v 80 -o reads/$file_name");

# Assemble
chdir("assemblies");
system("~/installations/VelvetOptimiser/VelvetOptimiser.pl --s 37 --e 81 --x 4 --t 1 -f '-shortPaired -fastq.gz -separate ../reads/$file_name\_100_250_1.fq.gz ../reads/$file_name\_100_250_2.fq.gz' -p $file_name");

# Clean large tmp files
my $best_assembly = `tail -1 $file_name\_logfile.txt`;
chomp $best_assembly;

rename "$best_assembly/contigs.fa", "$file_name\_contigs.fa";
remove_tree($best_assembly);

exit(0);
