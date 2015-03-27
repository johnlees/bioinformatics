#!/usr/bin/perl -w

package compare_variants;

use strict;
use warnings;

use Bio::SeqIO;

my $blastn_location = "blastn";
my $bcftools_location= "bcftools";

my %flip = (
   "A" => "T",
   "T" => "A",
   "G" => "C",
   "C" => "G");

sub classify_var($$$$)
{
   my (@variants) = @_;

   # Any site > 1 base must be an indel
   my ($type, $match);

   # Check for indels
   my @del_bases;
   for (my $i = 0; $i<=3; $i+=2)
   {
      if (length($variants[$i]) > 1 || length($variants[$i+1]) > 1)
      {
         $type = "INDEL";
         $del_bases[$i/2] = decompose_indel($variants[$i], $variants[$i+1]);
      }
      else
      {
         $del_bases[$i/2] = "";
      }

   }

   # Classify snps
   if (!defined($type))
   {
      $type = "SNP";

      my ($q_ref, $q_alt, $s_ref, $s_alt) = @variants;
      if ($q_ref eq $s_ref && $q_alt eq $s_alt)
      {
         $match = "EXACT";
      }
      elsif ($q_ref eq $s_alt && $q_alt eq $s_ref)
      {
         $match = "SWITCHED";
      }
      # This distinction can't really be made for A/T and G/C SNPs
      elsif (flip_strand($q_ref) eq $s_ref && flip_strand($q_alt) eq $s_alt)
      {
         $match = "EXACT_FLIP";
      }
      elsif (flip_strand($q_ref) eq $s_alt && flip_strand($q_alt) eq $s_ref)
      {
         $match = "SWITCH_FLIP";
      }
      else
      {
         $match = "MISMATCH";
      }
   }
   # Classify indels
   else
   {
      if ($del_bases[0] eq $del_bases[1])
      {
         $match = "EXACT";
      }
      elsif (flip_strand($del_bases[0]) eq $del_bases[1])
      {
         $match = "EXACT_FLIP";
      }
      elsif (flip_strand($del_bases[0], "reverse") eq $del_bases[1])
      {
         $match = "SWITCH_FLIP";
      }
      else
      {
         $match = "MISMATCH";
      }
   }

   # Work out match status. Assembly strands/contigs may be flipped!
   return($type, $match);
}

# Reverse complement a sequence
sub flip_strand($;$)
{
   my ($sequence, $direction) = @_;

   # Set whether to also reverse sequence
   if (!defined($direction))
   {
      $direction = "forward";
   }

   # Decompose string into bases
   my @bases = split("", $sequence);
   my @flipped;

   # Flip each base, and append to either start or end of array
   foreach my $base (@bases)
   {
      if ($direction eq "reverse")
      {
         unshift(@flipped, $flip{$base});
      }
      else
      {
         push(@flipped, $flip{$base});
      }
   }

   # Join bases back into a string and return
   my $flipped_sequence = join("", @flipped);
   return($flipped_sequence);
}

# Needs left aligned INDELS. Use bcftools norm
sub decompose_indel($$)
{
   my ($ref, $alt) = @_;

   # Treat everything as deletions i.e. move longer alts to be the ref
   if (length($alt) > length($ref))
   {
      my $tmp_ref = $ref;
      $ref = $alt;
      $alt = $tmp_ref;
   }

   my $del_bases;
   # INDEL could go at either end, as either end of ref sequence matches alt
   if (substr($ref, 0, length($alt)) eq substr($ref, length($ref) - length($alt), length($alt)))
   {
      # Homopolymer indels will hit this. They need to be treated
      # as below
      if (length($ref) == 2*length($alt))
      {
         $del_bases = substr($ref, length($alt));
      }
      else
      {
         $del_bases = substr($ref, length($alt), length($ref) - 2*length($alt));
      }
   }
   # For left aligned indels, can just take bases of ref past alt
   else
   {
      $del_bases = substr($ref, length($alt));
   }

   #DEBUG
   #print STDERR "$ref,$alt,$del_bases\n";
   return($del_bases);
}

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

# blastn of all v ref. Returns reference to hash of variant and new position
sub blastn_ref($$)
{
   my ($query_file, $ref_file) = @_;

   my %blast_scores;

   # blastn of all sequences in query with all in subject. Output is tab
   # sepearated query_id, subject_id, e_value and score for each pairwise
   # comparison
   my $blastn_command = "$blastn_location -subject $ref_file -query $query_file -outfmt \"6 qseqid sstart send evalue score\"";
   my $blast_result = `$blastn_command`;

   # Return data in a hash
   foreach my $pair (split("\n", $blast_result))
   {
      my ($q_id, $start, $end, $evalue) = split("\t", $pair);

      $blast_scores{$q_id}{start} = $start;
      $blast_scores{$q_id}{end} = $end;
      $blast_scores{$q_id}{evalue} = $evalue;
   }

   return(\%blast_scores);
}

# Extracts windows of a specified size around a specified list of variants,
# creates a multi-fasta file of them
sub variant_windows($$$$;$)
{
   my ($window_size, $variant_list, $sequence_file, $output_file, $before) = @_;

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
         my ($contig, $position, $ref, $alt, @samples) = split(",", $variant);

         if ($contig eq $sequence->display_id())
         {
            # Remove from variant list if found in this contig
            push (@found_variants, $variant);

            # Variant length. 1 for a SNP, > 1 for an insertion, < 1 for
            # a deletion
            my $var_length = length($alt) - length($ref) + 1;
            my $position_start = $position - 1;
            my $position_end;

            # SNPs
            if ($var_length == 1)
            {
               $position_end = $position + $var_length;
            }
            # Deletions
            elsif ($var_length < 1)
            {
               $position_end = $position - $var_length + 2;
            }
            # Insertions
            else
            {
               $position_end = $position + length($ref);
            }

            # Extract the window, then write it to the output
            my $window;
            if ($before)
            {
               $window = extract_bp_before($sequence, $window_size, $position_start);
            }
            else
            {
               $window = extract_bp_window($sequence, $window_size, $position_start, $position_end);
            }

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

   # Window size must be even, as we will extract a region symmetric around the
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

# Extracts sequence of a window size around a position, returning a Bio::Seq
# object of this
sub extract_bp_before($$$)
{
   my ($sequence, $seq_size, $position_start) = @_;

   my $context;
   if ($sequence->length() < $seq_size || $position_start < ($seq_size - 1) )
   {
      print STDERR "Cannot extract $seq_size bp before $position_start as sequence too small\n";
   }
   else
   {
      my $context_sequence = $sequence->subseq($position_start-$seq_size, $position_start);
      $context = Bio::Seq->new( -seq => $context_sequence);
   }

   return($context);
}

# Takes a list of paired vcfs and reference fastas and creates a multifasta of
# windows around the variants
sub create_blastn_input($$$;$$)
{
   my ($vcfs, $refs, $out_prefix, $filter, $type) = @_;

   my $i = 1;
   foreach my $vcf (@$vcfs)
   {
      my $variant_lists = extract_vcf_variants($vcf, $filter, $type);
      variant_windows(300, $variant_lists, shift(@$refs), "$out_prefix.$i.fa");

      $i++;
   }

}

# Get a list of variant sites from a vcf
sub extract_vcf_variants($;$$)
{
   my ($vcf_in, $filter, $type) = @_;

   # Tab delimited contig, position list of all variants in vcf
   my $bcftools_command;
   if (defined($type) && defined($filter))
   {
      $bcftools_command = "$bcftools_location view -f $filter -v $type $vcf_in | $bcftools_location query -f '%CHROM\t%POS\t%REF\t%ALT\n' -";
   }
   else
   {
      $bcftools_command = "$bcftools_location view $vcf_in | $bcftools_location query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' -";
   }
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
