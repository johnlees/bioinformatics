#!/usr/bin/perl -w

package quake_wrapper;

use strict;
use warnings;

use threads;
use File::Spec;
use Cwd;

my $quake_location = "/nfs/users/nfs_j/jl11/software/bin/quake/quake.py";

# Error corrects fastq files using quake
sub quake_error_correct($$$)
{
   my ($reads, $kmer_size, $threads) = @_;
   my (%corrected_reads, %quake_symlinks);

   mkdir "quake" || die("Could not create quake directory: $!\n");

   # Prepare a file with the locations of the fastq file pairs for input to
   # quake
   my $quake_input_file_name = "quake_reads.txt";
   open (QUAKE, ">quake/$quake_input_file_name") || die("Could not open $quake_input_file_name for writing: $!");

   # Quake outputs where the read files originally were, which isn't
   # necessarily writable. Create symlinks instead
   my $cwd = getcwd();
   foreach my $sample (keys %$reads)
   {
      foreach my $direction (keys %{$$reads{$sample}})
      {
         my $read_location = $$reads{$sample}{$direction};
         my ($volume ,$directories, $file) = File::Spec->splitpath($read_location);

         $quake_symlinks{$sample}{$direction} = $cwd . "/quake/$file";
         symlink $read_location, $quake_symlinks{$sample}{$direction};
      }
      print QUAKE join(" ", $quake_symlinks{$sample}{"forward"} , $quake_symlinks{$sample}{"backward"}) . "\n";
   }

   close QUAKE;

   # Run quake
   my $quake_command = "cd quake && $quake_location -f $quake_input_file_name -k $kmer_size -p $threads &> quake.log";
   system($quake_command);

   # Set paths of corrected reads
   foreach my $sample (keys %$reads)
   {
      foreach my $read_direction (keys %{$$reads{$sample}})
      {
         my ($volume ,$directories, $file) = File::Spec->splitpath($$reads{$sample}{$read_direction});
         $file =~ m/^(.+)\.(fastq|fq)$/;
         $corrected_reads{$sample}{$read_direction} = "$cwd/quake/$1.cor.$2";
      }
   }

   return(\%corrected_reads);
}

# Parses the read file names and sample names from the input read file.
# Returns reference to an array of sample names, and a hash of read
# locations.
sub parse_read_file($$)
{
   my ($read_file, $no_decompress) = @_;

   my (%read_locations, %decompress);
   my @sample_names;
   open(READS, $read_file) || die("Could not open $read_file: $!\n");

   while (my $read_pair = <READS>)
   {
      chomp($read_pair);
      my ($sample_name, $forward_read, $backward_read) = split("\t", $read_pair);
      push(@sample_names, $sample_name);

      # Decompress reads if needed, spawing a thread to do so
      if (defined($no_decompress) && !$no_decompress)
      {
         if ($forward_read =~ /\.gz$/)
         {
            $decompress{$sample_name}{"forward"} = threads->create(\&decompress_fastq, $forward_read);
         }
         if ($backward_read =~ /\.gz$/)
         {
            $decompress{$sample_name}{"backward"} = threads->create(\&decompress_fastq, $backward_read);
         }
      }

      # Store in hash of hashes
      $read_locations{$sample_name}{"forward"} = $forward_read;
      $read_locations{$sample_name}{"backward"} = $backward_read;
   }

   # Wait for all decompression threads to complete, and overwrite read
   # locations with new fastq paths
   foreach my $sample (keys %decompress)
   {
      foreach my $direction (keys %{$decompress{$sample}})
      {
         $read_locations{$sample}{$direction} = $decompress{$sample}{$direction}->join();
      }
   }

   close READS;
   return(\@sample_names, \%read_locations);
}

sub decompress_fastq($)
{
   my ($fastq) = @_;
   my $decompressed_location;
   my $cwd = getcwd();

   print STDERR "Decompressing $fastq\n";
   my ($volume ,$directories, $file) = File::Spec->splitpath($fastq);

   $file =~ m/^(.+\.)(fastq|fq)\.gz/;
   $decompressed_location = "$cwd/$1$2";
   system("gzip -d -c $fastq > $decompressed_location");

   return($decompressed_location);
}

1;

