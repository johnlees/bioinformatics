#!/usr/bin/perl -w

#****************************************************************************************#
#* resistance_extract.pl                                                                *#
#*                                                                                      *#
#****************************************************************************************#

use strict;

use Getopt::Long;

#* Gets input directory supplied by -i and reads the names of all contained files
my ($inputdir, $outputfile, $boolean);
GetOptions ("input|i=s"  => \$inputdir,
			"output|o=s" => \$outputfile,
			"boolean"    => \$boolean
		   ) or die($!);
		   
opendir (DIR, $inputdir) or die $!;		   
my @filelist = readdir(DIR);
closedir DIR;

open (OUTFILE, ">$outputfile");

# Specific genes to search for
my @boolean_array = ("catB9","dfrA1","strA","strB","sul1","sul2","dfrA18","floR");

# Header for boolean file
if ($boolean)
{
   foreach my $specific_gene (@boolean_array)
   {
      print OUTFILE ",$specific_gene";
   }
   print OUTFILE "\n";
}

foreach my $file (@filelist)
{
	# Look for resistance coverage summary files
	my $filepath = $inputdir . "/" . $file;
	if ($file =~ /^(\d+_\d+#\d+)_resistances_coverage\.txt$/)
	{
		my $sample_name = $1;
		
		open(STATSFILE, $filepath);
		# Do away with header row (gene, min, max, median, mean, std)
		my $linein = <STATSFILE>;
		
		my @genes;
		while ($linein = <STATSFILE>)
		{
			chomp($linein);
			my @stats_fields = split(/\s+/, $linein);
			
			# For each resistance gene the sample has, push the gene name into an array
			my $gene_name;
			if ($stats_fields[0] =~ /^(.+?)_/ && $stats_fields[3] > 5)
			{
			   $gene_name = $1;
			   push (@genes, $gene_name);
			}
			
		}
		
		close STATSFILE;
		
		print OUTFILE "$sample_name";
		
		# Print all genes
		if (!$boolean)
		{
		   foreach my $gene (@genes)
		   {
		      print OUTFILE ",$gene";
		   }
		}
		# Print only specified genes
		else
		{
		   # World's worst search algorithm
		   foreach my $specific_gene (@boolean_array)
		   {
		      my $i = 0;
		      my $found = 0;
		      while (!$found && $i < scalar(@genes))
		      {
		         if ($specific_gene eq $genes[$i])
		         {
		           $found = 1;
		         }
		         $i++;
		      }
		      
		      if ($found)
		      {
		         print OUTFILE ",+"
		      }
		      else
		      {
		         print OUTFILE ",-"
		      }
		   }
		}
		
		print OUTFILE "\n";
		
	}
}

close OUTFILE;

exit(0);