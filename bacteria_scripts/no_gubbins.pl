#!/usr/bin/perl -w

#****************************************************************************************#
#* no_gubbins.pl                                                                        *#
#* Takes the tab file produced by gubbins, parses it, and then changes all the regions  *#
#* identified as being due to recombination to Ns so they aren't considered during tree *#
#* construction. Uses Bio::AlignIO to read in and write to alignment files              *#
#****************************************************************************************#

use strict;

use Getopt::Long;

# BioPerl modules
use Bio::AlignIO;
use Bio::SeqIO;

use Data::Dumper;

my $help_message = <<END;
Usage: no_gubbins.pl -a <alignment_file> -t <tab_file> <options>
Removes the recombination regions deduced by gubbins, leaving the full alignment
with recombinant sites replaced by N, so they won't be used in the phylogeny

	-a, --align    The original alignment file gubbins was run on, in multi-
	               fasta format
	-t, --tab      The tab file output by gubbins, which lists the areas of
	               recombination to remove
	-o, --output   Name for the output multiple alignment, which is in multi-
	               fasta format. Defaults to no_gubbins.aln  
	-b, --begin    When only a region is required, the site this region starts
	               at
	-e, --end      The end of a region
	-h, --help     Displays this help message
END

#****************************************************************************************#
#* Subroutines                                                                          *#
#****************************************************************************************#

# Reads a tab file into a hash, sorts the arrays numerically for each hash key, and
# returns a reference to the hash
sub ReadGubbinsTab ($)
{
   my ($tab_file) = @_;
   
   my %tab_hash;

   open (TAB, $tab_file) or die("Couldn't open tab file $tab_file\n");
   
   # Read in tab data, which is in a fairly special format
   # Would be nice to use Text::CSV, but as far as I can tell this isn't really practical
   while (my $tab_line = <TAB>)
   {
      # Read an entry
      my $node_line = <TAB>;
      my $neg_log_lik_line = <TAB>;
      my $colour_line = <TAB>;
      my $taxa_line = <TAB>;
      my $SNP_line = <TAB>;
      
      # Parse range line
      my ($range_start, $range_end);
      chomp($tab_line);
      
      if ($tab_line =~ /^FT\s+misc_feature\s+(\d+)\.\.(\d+)$/)
      {
         $range_start = $1;
         $range_end = $2;
      }
      
      # Parse taxa line
      my @taxa_split;
      chomp($taxa_line);
      
      if ($taxa_line =~ /^FT\s+\/taxa="(.+)"$/)
      {
         @taxa_split = split(/\s+/, $1);
         
         if (scalar(@taxa_split) > 1)
         {
            shift(@taxa_split);
         }
         
      }
      
      # Put the range into the hash for each taxa specified
      foreach my $taxa (@taxa_split)
      {
         push(@{$tab_hash{$taxa}}, "$range_start,$range_end");
      }      
      
   }
   
   close TAB;
   
   # Sort by region start field
   foreach my $taxa (keys %tab_hash)
   {
      @{$tab_hash{$taxa}} = sort by_numeric_range @{$tab_hash{$taxa}};
   }
   
   # Rather than thinking hard about overlapping regions later, a better way might be
   # to collapse them all here so there are no overlapping regions specified in the array
   
   return(\%tab_hash);
   
}

# Sorts numerically ascending by the first field before the comma
sub by_numeric_range
{
	($a =~ /^(\d+),/)[0] <=> ($b =~ /^(\d+),/)[0];
}


#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* Gets input parameters
my ($alignment_file, $tab_file, $output_file, $region_start, $region_end, $help);
GetOptions ("align|a=s"  => \$alignment_file,
            "tab|t=s"    => \$tab_file,
            "output|o=s" => \$output_file,
            "begin|b=o"  => \$region_start,
            "end|e=o"    => \$region_end,
            "help|h"     => \$help
		   ) or die($help_message);
		   
# Check necessary files exist
if (defined($help))
{
   print $help_message;
}
elsif (!defined($alignment_file) || !defined($tab_file) || !-e $alignment_file 
                                                                         || !-e $tab_file)
{
	print ("Input alignment or tab file does not exist!\n\n");
	print $help_message;
}
else
{
   my $region_only = 0;
   
   # Check input params
   if ((defined($region_start) xor defined($region_end)))
   {
      print ("Both a valid start and end must be specified, or neither\n\n");
	  print $help_message;
   }
   elsif (defined($region_start) && ($region_start > $region_end))
   {
      print ("The start position of the region must be before the end value\n");
   }
   else
   {    
      # Make a default output file name
      if (!defined($output_file))
      {
         $output_file = "no_gubbins.aln";
      }
      
      print "Removing gubbins...\n\n";
      
      # Parse tab
      my $tab_data = ReadGubbinsTab($tab_file);
      
      # Good log output, uncomment if encountering problems
      #print Dumper($tab_data); 
      
      print("Tab file $tab_file read in successfully\n");
      print("Reading input alignment\n...\n");
   
      # Open alignment input and output. Alignment AlignI or SimpleAlign are objects in
      # AlignIO, which handles reading, writing and format only. Sequences are then
      # objects in AlignI or SimpleAlign objects
      my $alignment_in = Bio::AlignIO->new( -file   => "<$alignment_file",
                                            -format => "fasta" ) || die ($!);
      my $aln = $alignment_in->next_aln();
      
  	  my $alignment_out = Bio::AlignIO->new(-file   => ">$output_file",
                                            -format => "fasta",
                                            -flush  => 0 ) || die ($!);
      my $aln_out = new Bio::SimpleAlign();
                                            
      print("\nInput alignment read succesfully\n\n");                      
   
      # Use bioperl to do replacement serially
      foreach my $taxa (sort keys %$tab_data)
      {
         print("Removing recombination in alignment of $taxa\n");
         
         # Get array of sites to remove sorted by start value
         # For the given taxa's sequence
         foreach my $seq ( $aln->each_seq_with_id($taxa) )
         {
            my $i = 1;
            my $new_seq = '';
            
            while (1)
            {
               # Find start and end of next region, include sequence where not masked and
               # replace with the correct number of Ns otherwise
               my $next_range = shift(@$tab_data{$taxa});
               my ($start, $end);
               
               if (defined($next_range))
               {
                  my @start_and_end = split(',', $next_range);
                  $start = $start_and_end[0];
                  $end = $start_and_end[1];
               }
               else
               {
                  last;
               }
               
               # Bear in mind regions may be overlapping, so check for current pointer
               # position in case it's already in the region
               # e.g. regions overlap: 1-4 and 2-6 should mask out 1-6
               # e.g. region covered already: 1-4 and 2-3 should mask out 1-4
               if ($i < $start)
               {
                  $new_seq .= $seq->subseq($i, $start-1);
                  $new_seq .= "N"x($end-$start+1);
                  $i = $end + 1;
               }
               elsif ($end>$i)
               {
                  $new_seq .= "N"x($end-$i+1);
                  $i = $end + 1;
               }
               
            }
            
            # In case no sequence was replaced (e.g. in the reference)
            if ($new_seq eq '')
            {
               $new_seq = $seq->seq();
            }
            else
            {
               # complete end of sequence after final N
               $new_seq .= $seq->subseq($i, $seq->length());
            }
            
            # Make a sequence object, then add it to the simplealign object
            my $seq_write = Bio::LocatableSeq->new( -id => $taxa,
                                                    -seq => $new_seq);
            
            $aln_out->add_seq($seq_write);
            
         }
      }
      
      print ("Writing alignment out\n\n");
      
      # Print only the slice of alignment requested
      if (defined($region_start))
      {
         $alignment_out->write_aln($aln_out->slice($region_start, $region_end));
      }
      else
      {
         $alignment_out->write_aln($aln_out);
      }
      
      print ("Done!\n");
   }
   
}

exit(0);