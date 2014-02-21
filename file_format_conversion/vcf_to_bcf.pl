#!/usr/bin/perl -w

#****************************************************************************************#
#* vcf_to_bcf.pl                                                                        *#
#*                                                                                      *#
#****************************************************************************************#

use strict;

use Getopt::Long;

#* Gets input directory supplied by -i and reads the names of all contained files
my ($inputdir);
GetOptions ("input|i=s"  => \$inputdir
		   ) or die($!);
		   
opendir (DIR, $inputdir) or die $!;		   
my @filelist = readdir(DIR);
closedir DIR;

foreach my $file (@filelist)
{
   if ($file =~ /^(\d+_\d+#\d+)\.\d+\.mpileup\.unfilt\.vcf\.gz$/)
   {
      my $prefix = $1;
      my $command = "bcftools view -Sb -D chromosome.txt $inputdir/$file > $prefix.bcf";
      
      system($command);
   }
}