#!/usr/bin/perl -w
#
use strict;
use warnings;

use Getopt::Long;

use Bio::SeqIO;
use Bio::Tools::GFF;

#
# Globals
#
my $query_fasta_suffix = "hsdS_query.fa";

my $N_term_ref = "Ntermini.fasta";
my $C_term_ref = "Ctermini.fasta";

my $N_blast_suffix = "Ntermini.blast.out";
my $C_blast_suffix = "Ctermini.blast.out";
my $blast_err = "blast.err";

my $blast_evalue = "10e-6";

my %allele_map = ("Aa" => "A",
                  "Ab" => "B",
                  "Ac" => "E",
                  "Ba" => "D",
                  "Bb" => "C",
                  "Bc" => "F");

my $help_message = <<HELP;
Usage: ./ivr_typer.pl [--map|--assembly] <options>

Given S. pneumo read (via mapping mode) or annotated assembly (via assembly mode),
returns information on the likely allele type for the ivr/hsd R-M system locus

   Options
   --map                 Infer by mapping reads to reference alleles
   --assembly            Infer by annotation of gene, then a blast with reference
                         genes

   --ref_dir             Directory containing sequences of reference alleles
                         (default ./)

   <map mode>
   -f, --forward_reads   Location of forward reads (fastq(.gz))
   -r, --reverse_reads   Location of reverse reads (fastq(.gz))

   <assembly mode>
   -a, --annotation      Annotation of the -a assembly in gff v3.
   -b, --batch           Batch mode. Interprets the annotation option as a
                         list of gff files instead of a single file. A second
                         (tab separated) field with samples names can be used
                         STDOUT can be piped to get list of alleles by sample


   -h, --help            Displays this message

HELP

# Need to remove contigs with no features due to failings of Bio::Tools::GFF
# used with attach_seqs
# Should all be sorted by contig:position
sub pre_process_gff($$)
{
   my ($file_in, $file_out) = @_;

   my $fasta_tmp = "tmp.fa";

   if (-e $fasta_tmp)
   {
      die("$fasta_tmp already exists. Will not overwrite\n");
   }

   my $gff_in = Bio::Tools::GFF->new(-file => $file_in,
                                     -gff_version => 3) || die ("Could not open $file_in as gff v3: $!");

   my @contigs;
   my $last_contig = "";

   while (my $feature = $gff_in->next_feature())
   {
      if ($feature->seq_id() ne $last_contig)
      {
         $last_contig = $feature->seq_id();
         push(@contigs, $last_contig);
      }
   }

   my @contig_sequences = $gff_in->get_seqs();

   $gff_in->close();

   my $fasta_out = Bio::SeqIO->new( -file   => ">$fasta_tmp",
                                    -format => "fasta") || die("Could not write to $fasta_tmp: $!");

   foreach my $contig (reverse @contig_sequences)
   {
      foreach my $included_contig (@contigs)
      {
         if ($contig->display_id() eq $included_contig)
         {
            $fasta_out->write_seq($contig);
         }
      }
   }

   $fasta_out->close();

   #system("echo \"##FASTA\" | cat $gff_tmp - $fasta_tmp > $file_out");
   open(GFF, "$file_in") || die("Could not open $file_in\n");
   open(FASTA, "$fasta_tmp") || die("Could not open $fasta_tmp\n");
   open(OUT, ">$file_out") || die("Could not write to $file_out\n");

   while (my $gff_line = <GFF>)
   {
      if ($gff_line eq "##FASTA\n")
      {
         last;
      }
      else
      {
         print OUT $gff_line;
      }
   }
   print OUT "\n##FASTA\n";
   while (my $fasta_line = <FASTA>)
   {
      print OUT $fasta_line;
   }

   unlink $fasta_tmp;
   close OUT;

}

# Make the fasta file with the genes chosen by extract hsds. Include scores in
# sequence ids
sub make_query_fasta($$$$)
{
   my ($genes, $scores, $sequences, $f_out) = @_;

   my $sequence_out = Bio::SeqIO->new( -file   => ">$f_out",
                                       -format => "fasta") || die ($!);

   foreach my $sequence (@$sequences)
   {
      foreach my $feature ($sequence->get_SeqFeatures())
      {
         my @feature_ids = $feature->get_tag_values("ID");

         my $i = 0;
         foreach my $gene_id (@$genes)
         {
            my $score = $$scores[$i];

            my @hsds_ids = $gene_id->get_tag_values("ID");
            if ($hsds_ids[0] eq $feature_ids[0])
            {
               my $hsds_out = Bio::Seq->new( -seq => $feature->seq->seq(),
                                             -display_id => "$hsds_ids[0]_score$score");

               $sequence_out->write_seq($hsds_out);
            }

            $i++;
         }
      }
   }

   $sequence_out->close();
}

sub extract_hsds($)
{
   my ($annotation_file) = @_;

   my $tmp_annot_file = "annotation.tmp";
   pre_process_gff($annotation_file, $tmp_annot_file);

   my $gff_in = Bio::Tools::GFF->new(-file => $tmp_annot_file,
                                     -gff_version => 3) || die ("Could not open $tmp_annot_file as gff v3: $!");
   $gff_in->features_attached_to_seqs(1);

   # Hash of arrays, whose reference will be returned
   my %hsd_genes;
   my (@gene_order, @gene_strands);

   # Look through serially, find any hsd relevant genes
   while (my $feature = $gff_in->next_feature())
   {
      # Looking at products is more reliable than gene tags
      if ($feature->has_tag("product"))
      {
         my $strand = $feature->strand();
         push(@gene_strands, $strand);

         my $product = join(",", $feature->get_tag_values("product"));

         # This annotation is part of a type I R-M system
         if ($product =~ /type I restriction/i)
         {
            # Specificity subunit
            if ($product =~ /hsds/i || $product =~ /S protein/ || $product =~ /S subunit/ || $product =~ /chain S/)
            {
               push(@{ $hsd_genes{hsdS} }, $feature);
               push(@gene_order, "hsdS");
            }
            # Restriction subunit
            elsif ($product =~ /hsdr/i || $product =~ /R protein/ || $product =~ /R subunit/)
            {
               push(@{ $hsd_genes{hsdR} }, $feature);
               push(@gene_order, "hsdR");
            }
            # Methylation subunit
            elsif ($product =~ /hsdm/i || $product =~ /M protein/ || $product =~ /M subunit/)
            {
               push(@{ $hsd_genes{hsdM} }, $feature);
               push(@gene_order, "hsdM");

            }
         }
         # Possible creX/xerC gene, not too important to be specific here
         elsif($product =~ /recombinase/i)
         {
            push(@gene_order, "recombinase");
         }
         else
         {
            push(@gene_order, "-");
         }
      }
   }

   my @sequences = $gff_in->get_seqs();

   $gff_in->close();
   unlink $tmp_annot_file;

   # Now look through gene order, and try and find hsdS flanked by hsdM and
   # creX/xerC
   my (@hsds_genes, @scores);
   my $j = 0;

   for (my $i = 0; $i < scalar(@gene_order); $i++)
   {
      if ($gene_order[$i] eq "hsdS")
      {
         push(@hsds_genes, ${$hsd_genes{hsdS}}[$j]);

         # Score each hsdS gene. One point for each correctly annotated and
         # oriented gene in the locus. Score:1-6
         my $score = 1;

         if ($gene_order[$i-(2*$gene_strands[$i])] eq "hsdR" && $gene_strands[$i-(2*$gene_strands[$i])] == $gene_strands[$i])
         {
            $score++;
         }
         if ($gene_order[$i-(1*$gene_strands[$i])] eq "hsdM" && $gene_strands[$i-(1*$gene_strands[$i])] == $gene_strands[$i])
         {
            $score++;
         }
         if ($gene_order[$i+(1*$gene_strands[$i])] eq "recombinase" && $gene_strands[$i+(1*$gene_strands[$i])] == $gene_strands[$i])
         {
            $score++;
         }
         if ($gene_order[$i+(2*$gene_strands[$i])] eq "hsdS" && $gene_strands[$i+(2*$gene_strands[$i])] == $gene_strands[$i])
         {
            $score++;
         }
         if ($gene_order[$i+(3*$gene_strands[$i])] eq "hsdS" && $gene_strands[$i+(3*$gene_strands[$i])] == -$gene_strands[$i])
         {
            $score++;
         }

         push(@scores, $score);
         $j++;
      }
   }

   return(\@hsds_genes, \@scores, \@sequences);
}

# Run a blasn on subject and query
sub blastn($$$)
{
   my ($subject, $query, $output_file) = @_;

   my $blast_command = "blastn -subject $subject -query $query -evalue $blast_evalue -outfmt \"6 qseqid sallseqid evalue score\" > $output_file 2>> $blast_err";
   system($blast_command);

   # Parse output to extract best hits
   open(BLAST, "$output_file") || die("Could not open $output_file: $!\n");

   my $high_score = 0;
   my $best_gene_score = 0;
   my %blast_results;

   while (my $blast_line = <BLAST>)
   {
      chomp($blast_line);

      my ($query_id, $subject_id, $e_value, $score) = split("\t", $blast_line);

      #output highest scoring gene, not highest scoring blast hit
      $query_id =~ m/score(\d+)$/;
      my $gene_score = $1;

      if ($gene_score >= $best_gene_score)
      {
         $best_gene_score = $gene_score;

         if (!defined($blast_results{$query_id}{"score"}) || $score > $blast_results{$query_id}{"score"})
         {
            $blast_results{$query_id}{"score"} = $score;
            $blast_results{$query_id}{"hit"} = $subject_id;
         }
      }
   }

   close BLAST;

   return(\%blast_results);
}

# Pick the gene in both C and N terminals that gives the highest combined blast
# score
sub pick_hit_overlap($$)
{
   # Hash refs
   my ($hits1, $hits2) = @_;

   my %combined_scores;
   foreach my $query (keys %$hits1)
   {
      if (defined($$hits2{$query}{"score"}))
      {
         $combined_scores{$query} = $$hits1{$query}{"score"} + $$hits2{$query}{"score"};
      }
   }

   my ($hit1, $hit2, $query_id);
   my $high_score = 0;
   foreach my $query (keys %combined_scores)
   {
      if ($combined_scores{$query} > $high_score)
      {
         $hit1 = $$hits1{$query}{"hit"};
         $hit2 = $$hits2{$query}{"hit"};
         $query_id = $query;

         $high_score = $combined_scores{$query};
      }
   }

   return($query_id, $hit1, $hit2);
}

sub print_allele($$$$)
{
   my ($sample, $query, $N_term, $C_term) = @_;

   $N_term = seg_to_letter($N_term);
   $C_term = seg_to_letter($C_term);

   print join("\t", $sample, $query, "$N_term$C_term", $allele_map{"$N_term$C_term"}) . "\n";
}

sub seg_to_letter($)
{
   my ($term) = @_;

   if ($term =~ /^segment(.)$/)
   {
      $term = $1;
   }

   return($term);
}


#
# Main
#
my ($map, $assembly, $ref_dir, $forward_reads, $reverse_reads, $annotation_file_in, $batch, $help);
GetOptions ("map"       => \$map,
            "assembly"    => \$assembly,
            "ref_dir=s" => \$ref_dir,
            "forward_reads|f=s" => \$forward_reads,
            "reverse_reads|r=s" => \$reverse_reads,
            "annotation|a=s" => \$annotation_file_in,
            "batch|b" => \$batch,
            "help|h"     => \$help
		   ) or die($help_message);

if (!defined($ref_dir))
{
   $ref_dir = ".";
}

if (defined($help))
{
   print $help_message;
}
elsif (defined($map))
{

}
elsif (defined($assembly))
{
   if (!defined($annotation_file_in) || !-e $annotation_file_in)
   {
      die("Must set $annotation_file_in for assembly mode\n");
   }

   print STDERR "Using assembly\n";

   my (@annotation_files, %sample_names);
   if (defined($batch))
   {
      open(BATCH, $annotation_file_in) || die("Could not open $annotation_file_in\n");

      my $i = 0;
      while (my $batch_file_line = <BATCH>)
      {
         $i++;
         chomp $batch_file_line;

         my ($batch_file, $sample_name) = split("\t", $batch_file_line);

         if (-e $batch_file)
         {
            push(@annotation_files, $batch_file);
         }
         else
         {
            print STDERR "File $batch_file does not exist, skipping\n";
         }

         if (defined($sample_name) && $sample_name ne "")
         {
            $sample_names{$batch_file} = $sample_name;
         }
      }

      close BATCH;

      print STDERR "in batch mode for $i files\n";
   }
   else
   {
      # Single file
      push(@annotation_files, $annotation_file_in);
   }

   # Print header for output
   print STDERR "Most likely (highest scoring) allele per sample:\n";
   print join("\t", "sample", "query", "R6", "D39\n");

   foreach my $annotation_file (@annotation_files)
   {
      # This returns gff features, with attached sequence
      my ($hsds_genes, $scores, $sequences) = extract_hsds($annotation_file);

      # Process these objects into a multifasta to use as a blast query
      my ($query_fasta, $N_blast_out, $C_blast_out);
      if (defined($sample_names{$annotation_file}))
      {
         $query_fasta = $sample_names{$annotation_file} . ".$query_fasta_suffix";
         $N_blast_out = $sample_names{$annotation_file} . ".$N_blast_suffix";
         $C_blast_out = $sample_names{$annotation_file} . ".$C_blast_suffix";
      }
      else
      {
         $query_fasta = $query_fasta_suffix;
         $N_blast_out = $N_blast_suffix;
         $C_blast_out = $C_blast_suffix;
      }

      make_query_fasta($hsds_genes, $scores, $sequences, $query_fasta);

      # Run a protein blast for C and N termini
      my $N_term_blast = blastn("$ref_dir/$N_term_ref", $query_fasta, $N_blast_out);
      my $C_term_blast = blastn("$ref_dir/$C_term_ref", $query_fasta, $C_blast_out);

      my ($query, $N_hit, $C_hit) = pick_hit_overlap($N_term_blast, $C_term_blast);

      my $sample;
      if (defined($sample_names{$annotation_file}))
      {
         $sample = $sample_names{$annotation_file};
      }
      else
      {
         $sample = "-";
      }

      print_allele($sample, $query, $N_hit, $C_hit);
   }

}
else
{
   print STDERR "Mode must be set as one of map or assembly\n";
   print STDERR $help_message;
}

exit(0);

