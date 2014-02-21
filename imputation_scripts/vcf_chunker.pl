#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use Vcf;

sub vcfchunk ($$) {
	my ($input_vcf, $range) = @_;
	
	my ($chr, $from, $to);
	if ($range =~ /^chr(\d+)\:(\d+)\-(\d+)$/)
	{
		$chr = $1;
		$from = $2;
		$to = $3;
	}
	
	my $vcf_chunk = Vcf->new(file=>$input_vcf, region=>"$chr:$from-$to");
	
	return \$vcf_chunk;
}


# Set up environment
system("source ~/.bashrc");

my $help_message = "USAGE: --input --output --range chrX:100-200";

my ($input_file, $output_file, $range, $help);
GetOptions ("input=s"  => \$input_file,
            "output=s"    => \$output_file,
            "help|h"     => \$help
		   ) or die($help_message);

for (my $i = 1; $i < 14; $i++)
{
	open(VCFOUT, ">$i.$output_file.vcf") || die ("Couldn't open $i.$output_file.vcf for writing");
	
	my $start = ($i - 1)*5*1000000;
	my $end = $i*5*1000000;
	
	my $range = "chr20:$start-$end";
	
	my $vcf_chunk_ref = vcfchunk($input_file, $range);
	my $vcf_chunk = $$vcf_chunk_ref;
	
	$vcf_chunk->parse_header();
	print VCFOUT $vcf_chunk->format_header();
	
	my $prev;
    while (my $rec = $vcf_chunk->next_data_array())
    {
        if ( $$rec[1] < $start ) { next; }
        
        # skip duplicate positions
        if ( defined $prev && $prev eq $$rec[1] ) { next; }
        $prev = $$rec[1];

        print VCFOUT $vcf_chunk->format_line($rec);
    }
    $vcf_chunk->close();
    close(VCFOUT) || die ("Couldn't close $i.$output_file.vcf");
}

exit(0);