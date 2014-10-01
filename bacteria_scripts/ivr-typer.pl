#!/usr/bin/perl -w
#
use strict;
use warnings;

use Bio::SeqIO;

my $notes = <<NOTES;
Mapping:
./slice_fasta.pl
bwa index R6_spr0449.fa
bsub -o bwa.%J.o -e bwa.%J.e 'bwa mem ../R6_spr0449.fa 11861_4#46_1.fastq.gz 11861_4#46_2.fastq.gz > spr0449_mapping.sam'

bsub -q yesterday -R "select[mem>1000] rusage[mem=1000]" -M1000 -o sam_to_bam.log samtools sort -O bam -o spr0449_mapping.bam -T tmp spr0449_mapping.sam
bsub -q yesterday -R "select[mem>1000] rusage[mem=1000]" -M1000 -o sam_to_bam2.log samtools sort -n -O bam -o spr0449_mapping_qname.bam -T tmp spr0449_mapping.sam
rm spr0449_mapping.sam
samtools index spr0449_mapping.bam

samtools view -h -q 60 -F 0x904 -f 0x19 -o spr_5prime_downstream_unmapped.sam spr0449_mapping.bam
samtools view spr_5prime_downstream_unmapped.sam | cut -f 1 > spr_5prime_downstream_unmapped.qnames
OR
samtools view -h -q 60 -F 0x904 -f 0x19 spr0449_mapping.bam | cut -f 1 > spr_5prime_downstream_unmapped.qnames

may be best to then extract pairs from fastq_1 rather than qname sorted bam (difference is reverse complementing)?
wc -l spr_5prime_downstream_unmapped.qnames
80
zgrep -m 80 -F -A 1 -f spr_5prime_downstream_unmapped.qnames 11861_4\#46_1.fastq.gz
or from bam
samtools view spr0449_mapping.bam | grep -F -m 80 -f spr_5prime_downstream_unmapped.qnames | cut -f 10
(this is very fast on a qname sorted bam! The trick is maximum matches -m 80)
NOTES

open(READS, "downstream_reads.txt") || die("Could not open reads\n");

my $sequence_out = Bio::SeqIO->new( -file   => ">downstream_reads.fa",
                                    -format => "fasta") || die ($!);

my $i=1;
while (my $read = <READS>)
{
   chomp $read;

   my $new_read = Bio::Seq->new( -seq => $read,
                                 -display_id => "read_$i");
   $i++;

   $sequence_out->write_seq($new_read);
}

close READS;


exit(0);

