#!/usr/bin/perl -w

package assembly_pipeline;

use strict;
use warnings;

use File::Path qw(remove_tree);
use File::Copy;
use File::Spec;

use Bio::SeqIO;

# Program locations
my $spades_location = "/nfs/users/nfs_j/jl11/software/bin/spades.py";

# Program parameters
my $spades_kmers = "21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89";
my $spades_maxmem = 24; # in GiB

# Spades command, do not use directly
sub run_spades($$$$)
{
   my ($forward_reads, $reverse_reads, $output_dir, $threads) = @_;

   my $spades_command = "$spades_location -o $output_dir -1 $forward_reads -2 $reverse_reads --careful -t $threads -m $spades_maxmem -k $spades_kmers";

   system($spades_command);

   # Check output exists
   if (!-e "$output_dir/scaffolds.fasta" || !-e "$output_dir/contigs.fasta")
   {
      die("Assembly to $output_dir failed!\n");
   }

}

# Sub for spades assembly that works with tmp dirs
sub spades_assemble($$$$$)
{
   my ($forward_reads, $reverse_reads, $threads, $output_dir, $tmp_dir) = @_;

   if (!-d $tmp_dir)
   {
      mkdir $tmp_dir || die("Could not make $tmp_dir\n");
   }
   run_spades($forward_reads, $reverse_reads, $tmp_dir, $threads);

   copy("$tmp_dir/scaffolds.fasta", "$output_dir/scaffolds.fasta");
   copy("$tmp_dir/contigs.fasta", "$output_dir/contigs.fasta");

   remove_tree($tmp_dir);
}

# Filters spades output based on fasta headers
sub filter_contigs($$$$)
{
   my ($len_cutoff, $cov_cutoff, $input_file, $output_file) = @_;

   my $fasta_in = Bio::SeqIO->new(-file => $input_file,
                                  -format => 'fasta') || die($!);

   my $fasta_out = Bio::SeqIO->new(-file => ">$output_file",
                                  -format => 'fasta') || die($!);

   while (my $contig = $fasta_in->next_seq())
   {
      my ($length, $coverage);
      if ($contig->id =~ /^NODE_\d+_length_(\d+)_cov_(.+)/)
      {
         $length = $1;
         $coverage = $2;
      }

      if ($length >= $len_cutoff && $coverage >= $cov_cutoff)
      {
         $fasta_out->write_seq($contig);
      }
   }
}

# Don't call directly
sub run_improvement($$$$)
{
   my ($contigs_file, $forward_reads, $reverse_reads, $output_directory) = @_;

   # Get absolute paths
   $contigs_file = File::Spec->rel2abs($contigs_file);
   $forward_reads = File::Spec->rel2abs($forward_reads);
   $reverse_reads = File::Spec->rel2abs($reverse_reads);
   $output_directory = File::Spec->rel2abs($output_directory);

   # File names
   my @midway_improvements = ("scaffolds.filtered.fasta.scaffolded.filtered", "scaffolds.filtered.fasta.scaffolded.gapfilled.filtered", "scaffolds.scaffolded.gapfilled.length_filtered.fa");
   my $final_improvement = "scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa";
   my $renamed_improvement = "improved_assembly.fa";

   # Run improvement pipeline (requires pathogen softwarerc)
   if (!-d $output_directory)
   {
      mkdir $output_directory;
   }
   my $improve_command = "cd $output_directory && improve_assembly -a $contigs_file -f $forward_reads -r $reverse_reads -o $output_directory";
   system($improve_command);

   # Check file was successfully produced, otherwise die
   if (!-e "$output_directory/$final_improvement")
   {
      die("Improvement step in $output_directory failed!\n");
   }

   # Clear up files made part way through the improvement, and rename the
   # improved assembly to something more meaningful
   foreach my $midway_improvement (@midway_improvements)
   {
      unlink("$output_directory/$midway_improvement");
   }

   my $output_file = "$output_directory/$renamed_improvement";
   move("$output_directory/$final_improvement", $output_file);

   return($output_file);
}

# Use this for improvement. Fine to use scaffolds as input
sub improve_assembly($$$$$)
{
   my ($contigs_file, $forward_reads, $reverse_reads, $output_directory, $tmp_dir) = @_;

   if (!-d $tmp_dir)
   {
      mkdir $tmp_dir || die("Could not make $tmp_dir\n");
   }

   my $improved_assembly = run_improvement($contigs_file, $forward_reads, $reverse_reads, $tmp_dir);

   copy($improved_assembly, "$output_directory/improved_assembly.fa");

   remove_tree($tmp_dir);
}

# Unused, but still valid function (may be useful elsewhere?)
sub create_symlinks($$)
{
   my ($file_array, $link_directory) = @_;

   my @new_files;

   foreach my $file (@$file_array)
   {
      my ($volume,$directories,$file_name) = File::Spec->splitpath($file);
      my $new_name = "$link_directory/$file_name";

      symlink $file, $new_name;
      push(@new_files, $new_name);
   }

   return(\@new_files);
}

# The pipeline always outputs to a new directory 'annotation' in the pwd
sub run_annotation($$$$)
{
   my ($contigs_file, $sample_name, $genus, $threads) = @_;

   my $annotate_command = "annotate_bacteria -a $contigs_file --sample_name $sample_name --genus $genus --cpus $threads";
   system($annotate_command);

   if (-e "annotation/$sample_name.gff")
   {
      symlink "annotation/$sample_name.gff", "annotation.gff";
   }
   else
   {
      die("Annotation of $sample_name failed!\n");
   }
}

1;
