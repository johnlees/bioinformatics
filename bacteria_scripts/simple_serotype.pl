#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use mapping;

# Global settings
my $sort_mem = 2;
my $count_lim_default = 5;
my $num_serotypes_default = 3;
my $qual_cutoff_default = 20;
my $pileup_suffix = "pileup.txt";

# Reference location
my $references_default = "/lustre/scratch108/bacteria/jl11/reference_data/95_capsule_sequences.mfa";

my $usage_message = <<USAGE;
Usage: ./simple_serotype.pl -f forward.fastq -r reverse.fastq --id run_lane_tag <OPTIONS>

Maps input reads to reference serotypes and calculates maximally covered
sequences

   Options

   Required
   -f, --forward       Forward fasta reads
   -r, --reverse       Reverse fasta reads
   --id                Sample ID

   Optional
   --reference         Multi-fasta of capsule sequences
                       default: /lustre/scratch108/bacteria/jl11/reference_data/95_capsule_sequences.mfa
   -n, --num           Number of serotypes to display.
                       default: 3
   -c, --count         Number of reads mapping to each position to count read.
                       default: 5
   -q, --qual          Quality cutoff for mpileup
                       default: 20

   -h, --help          Shows this message

Requires: samtools, bcftools
At least one of: smalt, bwa, snap
USAGE

sub run_snap($$$$)
{
   my ($forward_reads, $reverse_reads, $sample_id, $reference_fasta) = @_;

   # Create reference index
   mapping::snap_index($reference_fasta);

   my $bam_name = "$sample_id.bam";
   my $snap_command = "$mapping::snap_location paired snap_index $forward_reads $reverse_reads -R '\@RG\\tID:$sample_id\\tSM:$sample_id'"
   . " -t 1 -= -so -sm $sort_mem -o $bam_name -om 1 -D 2 -F a &> snap_$sample_id.log";

   system($snap_command);

   return($bam_name);
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

# Read cmd line input
my ($forward_reads, $reverse_reads, $sample_id, $references, $num_serotypes, $count_lim, $qual_cutoff, $help);
GetOptions ("forward|f=s"  => \$forward_reads,
            "reverse|r=s"  => \$reverse_reads,
            "id=s"  => \$sample_id,
            "reference=s"  => \$references,
            "num|n=s"  => \$num_serotypes,
            "count|c=s"  => \$count_lim,
            "qual|q=s" => \$qual_cutoff,
            "help|h"     => \$help
		   ) or die($usage_message);

# Help?
if (defined($help))
{
   print $usage_message;
}
elsif (!-e($forward_reads) || !-e($reverse_reads))
{
   print STDERR "Could not find read files. Check they exist\n\n";
   print STDERR $usage_message;
}
elsif (!defined($sample_id))
{
   print STDERR "You must set the sample id\n\n";
   print STDERR $usage_message;
}
else
{
   # Set option defaults
   if (!defined($references))
   {
      $references = $references_default;
   }
   elsif (!-e $references)
   {
      print STDERR "Could not find references file $references";
      $references = $references_default;
   }

   if (!defined($num_serotypes))
   {
      $num_serotypes = $num_serotypes_default;
   }
   if (!defined($count_lim))
   {
      $count_lim = $count_lim_default;
   }
   if (!defined($qual_cutoff))
   {
      $qual_cutoff = $qual_cutoff_default;
   }

   # Map to reference sequences
   
   my $sero_bam = run_snap($forward_reads, $reverse_reads, $sample_id, $references);

   # Run and read pileup result. Count positions where at least a certain number of
   # reads map
   my $pileup_file = "$sample_id.$pileup_suffix";
   system("samtools mpileup -q $qual_cutoff -f $references $sero_bam > $pileup_file");
   open(PILEUP, $pileup_file) || die("Failed to open result $pileup_file\n");

   my %coverages;
   while (my $pileup_line = <PILEUP>)
   {
      chomp $pileup_line;
      my ($serotype, $pos, $allele, $count, @quals) = split("\t", $pileup_line);

      if ($count >= $count_lim)
      {
         $coverages{$serotype}++;
      }
   }

   close PILEUP;

   # Get reference sequence length from fasta index file, to normalise coverages
   if (!-e "$references.fai")
   {
      system("samtools faidx $references")
   }

   open(FAIDX, "$references.fai") || die("Could not open reference fasta index file $references.fai\nRun 'samtools faidx $references'\n");

   while (my $faidx_line = <FAIDX>)
   {
      chomp $faidx_line;

      my ($serotype, $length, @byte_pos) = split("\t", $faidx_line);
      if (defined($coverages{$serotype}))
      {
         $coverages{$serotype} /= $length;
      }
   }

   close FAIDX;

   # Sort to get the most covered serotypes
   my @top_serotypes = sort { $coverages{$b} <=> $coverages{$a} } keys(%coverages);
   print join("\t", splice(@top_serotypes, 0, $num_serotypes)) . "\n";

   my @top_coverages = @coverages{@top_serotypes};
   print join("\t", splice(@top_coverages, 0, $num_serotypes)) . "\n";

   # Clean up
   unlink($sero_bam, "$sero_bam.bai", $pileup_file,"snap_$sample_id.log");
}

exit(0);

