These are scripts I wrote (mostly in perl) to do some basic bioinformatics munging. They are no longer supported, but you are welcome to try using them. Archived readme follows

# README

To get a copy of them in your directory, run:

```git clone https://github.com/johnlees/bioinformatics/```

## assembly_scripts

Scripts generally to create and deal with de novo assemblies, or assembly only methods of calling

### annotate_vcf_from_gff.pl

Annotates a vcf with the features they fall in as well as gene names (using a gff file)

```./annotate_vcf_from_gff.pl <gff_in> <vcf_file>```

### mutate_fasta.pl
Simulates SNPs in a fasta file using the JC69 rate matrix and the specified mutation rate

```
./mutate_fasta.pl
Usage: ./mutate_fast.pl -f <fasta_file> -r <mutation_rate> -o <output_file> (-s <seed>) > mutations.txt

Creates SNPs in a multifasta using JC69 rate matrix
Options:

-f --fasta    Multi-fasta file of sequence to create SNPs in
-r --rate     Per-base mutation rate
-s --seed     Optional seed for random number generator
-o --output   Output file name

-h --help     This help message

Prints mutations tab separated to stdout and logs to stderr
```

### run_quake.pl

Simple wrapper to run quake on pairs of read files (fastq(.gz)) - which error corrects reads based on kmer frequencies. For the kmer to use see the help
This should be bsubbed with `bsub -o quake.%J.o -n4 -R "span[hosts=1]" -R "select[mem>3000] rusage[mem=3000]" -M3000 ./run_quake.pl <options>`

```
./run_quake.pl

One of the options was not defined. All are required
Usage: ./run_quake.pl -r <reads_file> -k <kmer_size> -t <threads>

Runs quake error correction on paired end reads provided as fastq(.gz) files.
Options:
-r --reads    Tab delimited file of fastq or fastq.gz file locations.
One line per sample. First column sample name, second
column forward reads, third column reverse reads.
It is best to give absolute paths
-k --kmer     The kmer size to use with quake. Set to
ceil[log(200*G)/log(4)] where G is the genome length
-t --threads  The number of threads to allow quake to use
-s --separate Error correct each sample separately. Probably required
for simulated data

-h --help     This help message

Output will be in cwd/quake
```

### reference_free_variant_caller.pl

Calls variants (at the moment only SNPs, but soon INDELs too) between two samples without needing to provide a reference. The output is an annotated vcf with just the variants

You'll need to provide the assembly of one of the samples as input (it is only used to place the variants, not in calling). See the extended help with -h for advice on how to prepare this.

bsub commands can be found in the extended help by running with the flag -h. I'd advise using 2 or 4 cores

You can use either cortex or sga to do this. Early tests suggests cortex has better power. Reads are error corrected with quake first

```
/reference_free_variant_caller.pl
Usage: ./reference_free_variant_caller.pl [--sga|--cortex] -a <assembly.fasta> -g <annotation.gff> -r <reads_file>

Uses cortex_var or sga to call variants between two samples without providing a reference
sequence

Options
--cortex            Use cortex for variant calling
--sga               Use sga for variant calling

-a, --assembly      Contigs of a de-novo assembly of one of the samples in
multi-fasta format. See help for more info.
-g, --annotation    Annotation of the -a assembly in gff v3.
-r, --reads         Tab delimited file of fastq or fastq.gz file locations.
One line per sample. First column sample name, second
column forward reads, third column reverse reads.
It is best to give absolute paths

-o, --output        Prefix for output vcfs

-c, --med-cov       Median coverage of reads over the assembly. Required by
sga, unused by cortex

--no-error-correct  Don't run quake on reads. Input reads are fed directly
into cortex
--separate-correct  Error correct each sample separately. Probably required
for simulated data

--dirty             Don't clean up temporary files

-h, --help          Shows more detailed help.

For cortex, requires: cortex_var, bcftools, quake and jellyfish.
For sga, requires: sga
See help for more details
```

## assembly_pipeline

Run `./batch_assemble.pl`

This takes a list of run_lane#tag and runs (over LSF):

* SPAdes assembler (using --careful mode which include hammer read error correction)
* improvement pipeline
* annotation pipeline

This script does not need to be bsub'd. The last three options can be used to continue failed jobs (e.g. with more memory)

```
./batch_assemble.pl
Usage: ./batch_assemble.pl <options> --lanes lane_file.txt

Submit assembly and annotation jobs across lsf
Takes a file which contains ids in the form run_lane#tag, one per line

Options
--lanes       A list, one per line, of sequencing ids of the form
run_lane#tag
--genus       The genus of the bacteria (used for annotation)

--output      Output directory (default ./)
--tmp         Directory to write temporary files to (default /tmp)
--threads     Number of threads to use (default 1)
--mem         Comma separated list of memory to submit each step with
Default: 4000,2000,1200

--nosymlinks  Do not create subdirectories with symlinks, as they
already exist
--improve     Run from improvement step
--annotate    Run from annotation step

--help      Displays this help message
```

## bacteria_scripts

Anything else

### compare_variants.pl
Does a blastn of nucleotides in a 300bp window around a variant. This can be used to confirm variants are the same when they have been called after mapping onto different reference sequences

```
./compare_variants.pl --vcf1 1.vcf.gz --vcf2 2.vcf.gz --ref1 ref1.fa --ref2 ref2.fa
```

### ivr-typer.pl

*Strep. pneumo. specific*

Extracts the allele type for the hsdS gene in the type I R-M system in S. pneumo which is one of six types (see Manso et al, Nat. Comms. 2014).

You can try and do this from assemblies, or using mapping. They work in different ways, but I think mapping is better

Run with: 
```
./ivr-typer.pl --map --ref_dir /lustre/scratch108/bacteria/jl11/pairs_analysis/RM-locus/mapping -f forward_reads.fastq.gz -r reverse_reads.fastq.gz > ivr_allele.txt
```

```
./ivr-typer.pl [--map|--assembly] <options>

Given S. pneumo reads (via mapping mode) or annotated assembly (via assembly mode),
returns information on the likely allele type for the ivr/hsd R-M system locus

Using mapping mode will take around 1Gb memory for bam sorting, and around 5 mins
of CPU time if mapping is required

STDOUT can be piped to get tab separated output with a header

Options
--map                 Infer by mapping reads to reference alleles
--assembly            Infer by annotation of gene, then a blast with reference
genes

--ref_dir             Directory containing sequences of reference alleles
(default ./)

<map mode>
-m, --mapping         A bam file which already has reads mapped to the D39 ref
OR (to create the bam)
-f, --forward_reads   Location of forward reads (fastq(.gz))
-r, --reverse_reads   Location of reverse reads (fastq(.gz))

<assembly mode>
-a, --annotation      Annotation of the -a assembly in gff v3.
-b, --batch           Batch mode. Interprets the annotation option as a
list of gff files instead of a single file. A second
(tab separated) field with samples names can be used
STDOUT can be piped to get list of alleles by sample


-h, --help            Displays this message
```

To make sense of the output of this script run on many samples you can use `ivr_format.pl`

### map_snp_call.pl
Very similar to `~sh16/scripts/multiple_mappings_to_bam.py`, developed originally to map two very similar sequences to a draft assembly. More features have been added recently

Designed for two samples at a time (c.f. reference_free_variant_caller.pl), but put in as many/few as you want.
Creates one large merged bam and calls from this. I doubt this will work with hundreds/thousands of samples

Advantages are:

* Uses new bcftools caller which should be better
* Calls are haploid not diploid
* Allows use of a more accurate prior for bacteria rather than the human default
* VCFs are annotated if they fall in a region with a feature. Frameshifting indels are annotated
* Can use snap for mapping, which is much quicker than bwa mem or smalt
* Works with multiple samples concurrently

Disadvantages:

* No pseudo-genome
* Not great scaling (use the same number of threads as samples, but would be better as job arrays). This should be ok if you write a script to bsub each sample
* Memory scaling of snap seems capricious
* Multiple cores is dodgy with snap

Future improvements:

* Allow use of GATK genotyper (quicker)
* Try and call large deletions

```
Usage: ./map_snp_call.pl -r <reference.fasta> -a <reference_annotation.gff> -r <reads_file> -o <output_prefix>

Maps input reads to reference, and calls snps

Options
-a, --assembly      Fasta file of reference to map to.
-g, --annotation    Annotation of the reference to annotate variants with.
-r, --reads         Tab delimited file of fastq or fastq.gz file locations.
One line per sample. First column sample name, second
column forward reads, third column reverse reads.
It is best to give absolute paths
-o, --output        Prefix for output files

-m, --mapper        Choose snap, bwa or smalt. Default snap
--gatk              Use gatk to improve bams (realign around indels)
(recommended!)

-t, --threads       Number of threads to use. Default 1
-p, --prior         Prior for number of snps expected. Default 1e-3

--dirty             Don't remove temporary files

-h, --help          Shows more detailed help.

Requires: samtools, bcftools
At least one of: smalt, bwa, snap
Optional: picard, gatk
```


