#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($inputdir);
GetOptions ("input|i=s"  => \$inputdir
		   ) or die($!);
		   
if (defined($inputdir))
{
   chdir($inputdir);
}

my @file_list = glob("*.plot");

foreach my $file (@file_list)
{
   if ($file =~ /^coverage\.(\d+_\d+#\d+)\..+\.bam\.plot/)
   {
      my $newfile = "$1.plot";
      system("cp $file $newfile");
   }
}

exit(0);
