#!/usr/bin/perl -w

#****************************************************************************************#
#* map_resistance.pl                                                                    *#
#* Runs map_resistome.py on all assemblies in the passed directory. Jobs are run        *#
#* serially and not in parallel to avoid causing reams of huge temporary files taking   *#
#* up all the disk space                                                                *#
#* Most things are hard coded...                                                        *#
#****************************************************************************************#

use strict;
use Getopt::Long;

#* 
#* Main
#*

#* Gets input directory supplied by -i and reads the names of all contained files
my $inputdir;
GetOptions ("input|i=s"  => \$inputdir,
		   ) or die($!);
		   
opendir (DIR, $inputdir) or die $!;		   
my @filelist = readdir(DIR);
closedir DIR;

my $prev_job_num;

foreach my $file (@filelist)
{
	# For each assembly, get the run and lane number
	if ($file =~ /^(\d+_\d+#\d+)_automated_velvet\.fa$/)
	{
		my $command;
		my $file_ref = $1;
		
		# Take the output from the previous submitted command, and bsub using the
		# condition that the previous job number must have completed before leaving the
		# PEND state
		if (defined($prev_job_num) && $prev_job_num =~ /<(\d+)>/)
		{
			$command = "bsub -w 'done($1)' -o $file_ref" . "_resistance.log ~sh16/scripts/map_resistome.py -c $inputdir/$file_ref" . "_automated_velvet.fa -f $inputdir/$file_ref" . "_1.fastq.gz -r $inputdir/$file_ref" . "_2.fastq.gz -g Curated_Database_20131022.fasta -o $file_ref" . "_resistances";
		}
		else
		{
			$command = "bsub -o $file_ref" . "_resistance.log ~sh16/scripts/map_resistome.py -c $inputdir/$file_ref" . "_automated_velvet.fa -f $inputdir/$file_ref" . "_1.fastq.gz -r $inputdir/$file_ref" . "_2.fastq.gz -g Curated_Database_20131022.fasta -o $file_ref" . "_resistances";
		}
		
		$prev_job_num = `$command`;
		print "$command\n";
		print "$prev_job_num\n\n";
	}
}

exit(0);

