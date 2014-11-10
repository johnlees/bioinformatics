#!/usr/bin/perl -w

package compare_variants;

use strict;
use warnings;

use Bio::SeqIO;

my $blastn_location = "blastn";
my $bcftools_location= "bcftools";

# list diff, union, intersection
sub list_compare($$$)
{
   my ($mode, $list1, $list2) = @_;

   my %count;
   my (@union, @intersection, @diff);

   foreach my $element (@$list1, @$list2)
   {
      $count{$element}++;
   }
   foreach my $element (keys %count)
   {
      push(@union, $element);

      if ($count{$element} > 1)
      {
         push(@intersection, $element);
      }
      else
      {
         push(@diff, $element);
      }
   }

   my $return_ref;
   if ($mode eq "diff")
   {
      $return_ref = \@diff;
   }
   elsif ($mode eq "intersection")
   {
      $return_ref = \@intersection;
   }
   elsif ($mode eq "union")
   {
      $return_ref = \@union;
   }
   else
   {
      die("Invalid list mode $mode\n");
   }

   return($return_ref);
}


# blastn of all v all. Returns reference to hash of sample by sample scores
sub blastn_pairwise($$)
{
   my ($query_file, $subject_file) = @_;

   my %blast_scores;

   # blastn of all sequences in query with all in subject. Output is tab
   # sepearated query_id, subject_id, e_value and score for each pairwise
   # comparison
   my $blastn_command = "$blastn_location -subject $subject_file -query $query_file -outfmt \"6 qseqid sallseqid evalue score\"";
   my $blast_result = `$blastn_command`;

   # Return data in a hash
   foreach my $pair (split("\n", $blast_result))
   {
      my ($q_id, $s_id, $e_val, $score) = split("\t", $pair);

      $blast_scores{$q_id}{$s_id} = $score;
   }

   return(\%blast_scores);
}

# Extracts windows of a specified size around a specified list of variants,
# creates a multi-fasta file of them
sub variant_windows($$$$)
{
   my ($window_size, $variant_list, $sequence_file, $output_file) = @_;

   # Input and output multifasta
   my $sequence_in = Bio::SeqIO->new( -file   => "<$sequence_file",
                                      -format => "fasta" ) || die ($!);
   my $sequence_out = Bio::SeqIO->new( -file   => ">$output_file",
                                       -format => "fasta") || die ($!);

   # Each fasta sequence, check all remaining variants
   my @found_variants;
   while (my $sequence = $sequence_in->next_seq())
   {
      foreach my $variant (@$variant_list)
      {
         my ($contig, $position, $ref, $alt) = split(",", $variant);

         if ($contig eq $sequence->display_id())
         {
            # Remove from variant list if found in this contig
            push (@found_variants, $variant);

            # Variant length. 1 for a SNP, > 1 for an insertion, < 1 for
            # a deletion
            my $var_length = length($alt) - length($ref) + 1;
            my ($position_start, $position_end);

            # SNPs
            if ($var_length == 1)
            {
               $position_start = $position - 1;
               $position_end = $position + $var_length;
            }
            # Deletions
            elsif ($var_length < 1)
            {
               $position_start = $position;
               $position_end = $position - $var_length + 2;
            }
            # Insertions
            else
            {
               $position_start = $position;
               $position_end = $position + 1;
            }

            # Extract the window, then write it to the output
            my $window = extract_bp_window($sequence, $window_size, $position_start, $position_end);

            $window->display_id($variant);
            $sequence_out->write_seq($window);
         }
      }
   }

   # Print a list of missed variants
   my $missed_variants = list_compare("diff", $variant_list, \@found_variants);
   foreach my $variant (@$missed_variants)
   {
      print STDERR "Variant $variant not found in input $sequence_file\n";
   }
}

# Extracts sequence of a window size around a position, returning a Bio:Seq
# object of this
sub extract_bp_window($$$;$)
{
   my ($sequence, $window_size, $position_start, $position_end) = @_;

   # Window size must be evn, as we will extract a region symmetric around the
   # variant not including it
   if ($window_size % 2 != 0)
   {
      $window_size--;
   }
   my $half_size = ($window_size)/2;

   my ($window, $start, $end);

   if ($sequence->length() < ($window_size + $position_end - $position_start))
   {
      print STDERR "Cannot extract $window_size bp around $position_start as sequence too small\n";
   }
   else
   {
      if (($position_end + $half_size) > $sequence->length())
      {
         # Extract from end backwards
         $start = $sequence->length() - ($window_size + $position_end - $position_start);
         $end = $sequence->length();
      }
      elsif (($position_start - $half_size) <= 0)
      {
         # Extract from start forwards
         $start = 1;
         $end = $window_size + ($position_end - $position_start);
      }
      else
      {
         $start = $position_start - $half_size + 1;
         $end = $position_end + $half_size - 1;
      }

      my $window_sequence = $sequence->subseq($start, $position_start) . $sequence->subseq($position_end, $end);
      $window = Bio::Seq->new( -seq => $window_sequence);
   }

   return($window);
}

# Takes a list of paired vcfs and reference fastas and creates a multifasta of
# windows around the variants
sub create_blastn_input($$$)
{
   my ($vcfs, $refs, $out_prefix) = @_;

   my $i = 1;
   foreach my $vcf (@$vcfs)
   {
      my $variant_lists = extract_vcf_variants($vcf);
      variant_windows(300, $variant_lists, shift(@$refs), "$out_prefix.$i.fa");

      $i++;
   }

}

# Get a list of variant sites from a vcf
sub extract_vcf_variants($)
{
   my ($vcf_in) = @_;

   # Tab delimited contig, position list of all variants in vcf
   my $bcftools_command = "$bcftools_location query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf_in";
   my $bcf_return = `$bcftools_command`;

   # Create a hash of variant details
   my @variant_list;
   foreach my $variant (split("\n", $bcf_return))
   {
      $variant =~ s/\t/,/g;
      push (@variant_list, $variant);
   }

   return(\@variant_list);
}

1;
