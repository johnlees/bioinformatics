#!/usr/bin/perl -w

package assembly_pipeline;

use strict;
use warnings;

use File::Path qw(remove_tree);
use File::Copy;
use File::Spec;
use Cwd;

use Bio::SeqIO;

# Program locations
my $spades_location = "/nfs/users/nfs_j/jl11/software/bin/spades.py";

# Program parameters
my $spades_kmers = "21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89";
my $spades_maxmem = 8; # in GiB

sub run_spades($$$$)
{
   my ($forward_reads, $reverse_reads, $output_dir, $threads) = @_;

   my $spades_command = "$spades_location -o $output_dir -1 $forward_reads -2 $reverse_reads --careful -t $threads -m $spades_maxmem -k $spades_kmers";

   system($spades_command);
}

sub spades_assemble($$$$)
{
   my ($forward_reads, $reverse_reads, $threads, $output_dir) = @_;

   my $tmp_dir = "/tmp/spades";

   mkdir $tmp_dir || die("Could not make $tmp_dir\n");
   run_spades($forward_reads, $reverse_reads, $tmp_dir, $threads);

   copy("$tmp_dir/scaffolds.fasta", "$output_dir/scaffolds.fasta");
   copy("$tmp_dir/contigs.fasta", "$output_dir/contigs.fasta");

   remove_tree($tmp_dir);
}

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
      if ($contig->id =~ /^NODE_\d+_length_(\d+)_cov_(.+)_ID_\d+$/)
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

# Fine to use scaffolds as input
sub run_improvement($$$$)
{
   my ($contigs_file, $forward_reads, $reverse_reads, $output_directory) = @_;

   # File names
   my @midway_improvements = ("scaffolds.filtered.fasta.scaffolded.filtered", "scaffolds.filtered.fasta.scaffolded.gapfilled.filtered", "scaffolds.scaffolded.gapfilled.length_filtered.fa");
   my $final_improvement = "scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa";
   my $renamed_improvement = "improved_assembly.fa";

   # Run improvement pipeline (requires pathogen softwarerc)
   my $improve_command = "improve_assembly -a $contigs_file -f $forward_reads -r $reverse_reads -o $output_directory";
   system($improve_command);

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

sub improve_assembly($$$$)
{
   my ($contigs_file, $forward_reads, $reverse_reads, $output_directory) = @_;

   my $tmp_dir = "/tmp/improvement";
   mkdir $tmp_dir || die("Could not make $tmp_dir\n");

   my $new_files = create_symlinks(\($contigs_file, $forward_reads, $reverse_reads), $tmp_dir);
   my $old_wd = getcwd();

   chdir($tmp_dir);
   my $improved_assembly = run_improvement(shift(@$new_files), shift(@$new_files), shift(@$new_files), $tmp_dir);
   chdir($old_wd);

   copy($improved_assembly, "$output_directory/improved_assembly.fa");

   remove_tree($tmp_dir);

}

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

   symlink "annotation/$sample_name.gff", "annotation.gff";
}

1;
