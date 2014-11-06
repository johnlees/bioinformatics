#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Based on Martin Kircher's script apply_parameters.py downloaded from http://cadd.gs.washington.edu/simulator
Used in:
"A general framework for estimating the relative pathogenicity of human genetic variants" by Kircher & Witten et al. 2014

Repurposed for use on 2Mb bacterial genomes by John Lees:
Parameters are estimated from progessiveCactus alignments rather than Ensembl Compara EPO pipeline
Chromosome options removed
Lower default number of events
No use of regions
Use GATK to produce fasta files from the vcfs

"""

import os,sys
from collections import defaultdict
from optparse import OptionParser
import traceback
import random
import mmap
import subprocess

java_location = "/software/bin/java"
gatk_location = "~/software/bin/GenomeAnalysisTK.jar"

# READ FASTA INDEX TO MEMORY
def read_fasta_index(filename):
  infile = open(filename)
  res = {}
  for line in infile:
    fields = line.split()
    if len(fields) == 5:
      cname,length,start,line,cline = fields
      res=cname,int(length),int(start),int(line),int(cline)
    else:
      sys.stderr.write('Error: Unexpected line in fasta index file: %s\n'%(line.strip()))
      sys.exit()
  infile.close()
  return res

parser = OptionParser()
parser.add_option("-n","--number_events", dest="number", help="Approximate number of substitution events to simulate (default 100)",type="int",default=100)
parser.add_option("-l","--log_file",dest="log_file",help="The location of the log file produced by hal2epo.pl",default="bacteria.log")
parser.add_option("-i", "--infile", dest="infile", help="Infile name of genome (prefix for fasta index, def reference_genome.fa)",default="reference_genome.fa")
parser.add_option("-o", "--outfile", dest="outfile", help="Prefix for output file (default STDOUT)",default=None)
parser.add_option("-f", "--fasta", action="store_true", dest="fasta_out", help="Produce a fasta file with the mutations (default false)",default=False)
parser.add_option("-s", "--seed", dest="seed", help="Seed for the random number generator (to generate the same mutations)",default=None)
(options, args) = parser.parse_args()

# OPEN OUTPUT FILE IF GIVEN, OTHERWISE PRINT ON SCREEN
if options.outfile != None:
  output = open(options.outfile+".vcf",'w')
else:
  output = sys.stdout

if options.fasta_out == True:
  if options.outfile == None:
    sys.stderr.write("To write a fasta out, an outfile -o must be specified")
  else:
    fasta_out_path = options.outfile+".fa"

if not os.path.exists(options.infile) or not os.path.exists(options.infile+".fai"):
  sys.stderr.write("Invalid path for genome and genome fasta index file.\n")
  sys.exit()

# Set seed if necessary
if options.seed != None:
  random.seed(options.seed)

fastaindex = read_fasta_index(options.infile+".fai")

if not os.path.exists(options.log_file):
  sys.stderr.write("Invalid path for log file.\n")
  sys.exit()

totalrefA = 0
totalrefC = 0
totalrefG = 0
totalrefT = 0

mut,total = 0,0
nAC,nAG,nAT,nCA,nCG,nCT,nGA,nGC,nGT,nTA,nTC,nTG = 0,0,0,0,0,0,0,0,0,0,0,0
totalCpG,mutCpG = 0,0
mCA,mCG,mCT,mGA,mGC,mGT = 0,0,0,0,0,0

insertionsizes,deletionsizes = defaultdict(int),defaultdict(int)

intervals = defaultdict(list)
interval_vals = {}

try:
  infile = open(options.log_file)
  line = infile.readline()
  while line != "":
    if line == "#A\tC\tG\tT\tCpGs\n":
      line = infile.readline()
      A,C,G,T,cCpG = map(int,line.split())
      totalrefA += A
      totalrefC += C
      totalrefG += G
      totalrefT += T
      totalCpG += cCpG
    elif line == "#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n":
      line = infile.readline()
      cmut,ctotal,cAC,cAG,cAT,cCA,cCG,cCT,cGA,cGC,cGT,cTA,cTC,cTG = map(int,line.split())
      mut+=cmut
      total+=ctotal
      nAC+=cAC
      nAG+=cAG
      nAT+=cAT
      nCA+=cCA
      nCG+=cCG
      nCT+=cCT
      nGA+=cGA
      nGC+=cGC
      nGT+=cGT
      nTA+=cTA
      nTC+=cTC
      nTG+=cTG
    elif line == "#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n":
      line = infile.readline()
      cmutCpG,ctotalCpG,cCA,cCG,cCT,cGA,cGC,cGT = map(int,line.split())
      mutCpG+=cmutCpG
      totalCpG+=ctotalCpG
      mCA+=cCA
      mCG+=cCG
      mCT+=cCT
      mGA+=cGA
      mGC+=cGC
      mGT+=cGT
    elif line == "##INSERTIONS\n":
      line = infile.readline() # HEADER
      line = infile.readline()
      while line != "" and not line.startswith('#'):
        fields = line.split()
        if len(fields) == 2:
          length,count=map(int,fields[:2])
          insertionsizes[length]+=count
        line = infile.readline()
    elif line == "##DELETIONS\n":
      line = infile.readline() # HEADER
      line = infile.readline()
      while line != "" and not line.startswith('#'):
        fields = line.split()
        if len(fields) == 2:
          length,count=map(int,fields[:2])
          deletionsizes[length]+=count
        line = infile.readline()
    else:
      line = infile.readline()
  infile.close()
except:
  exc_type, exc_value, exc_traceback = sys.exc_info()
  sys.stderr.write("%s\n"%str(exc_value))
  traceback.print_tb(exc_traceback)
  sys.stderr.write('Script terminated early. Printing current values.\n')

sum_obs = float(totalrefA+totalrefC+totalrefG+totalrefT)
total = float(total)
totalCpG = float(totalCpG)
sys.stderr.write("Total aligned positions: %d\n"%(sum_obs))
sys.stderr.write("Total mutations: %d/%d; Rate: %.8f/%.8f\n"%(mut,mutCpG,mut/total,mutCpG/totalCpG))
gmut = mut/total
gmutCpG = mutCpG/totalCpG

gmut_ = mut/sum_obs
gmutCpG_ = mutCpG/sum_obs
sys.stderr.write("Abs. mutation rate non-CpG: %.8f\tAbs. mutation rate CpGs: %.8f\n"%(gmut_,gmutCpG_))

dmut = (options.number*gmut_/(gmut_+gmutCpG_))/total
dmutCpG = (options.number*gmutCpG_/(gmut_+gmutCpG_))/totalCpG
dmutindel = options.number/sum_obs
sys.stderr.write("Likelihood altering a CpG: %.8f; Altering another base: %.8f\n"%(dmutCpG,dmut))

# RATES DETERMINED FROM BACTERIAL ANCESTOR USING PROGRESSIVECACTUS (CORRECTED FOR BASE COMPOSITION)

totalrefA = float(totalrefA)
totalrefC = float(totalrefC)
totalrefG = float(totalrefG)
totalrefT = float(totalrefT)

bfreqs = [totalrefA/total,totalrefC/total,totalrefG/total,totalrefT/total]

##SUBSTITUTION MATRIX NON-CpG
GTR_nonCpG = {'A':[             0 , nAC/totalrefA , nAG/totalrefA , nAT/totalrefA ],
              'C':[ nCA/totalrefC ,             0 , nCG/totalrefC , nCT/totalrefC ],
              'G':[ nGA/totalrefG , nGC/totalrefG ,             0 , nGT/totalrefG ],
              'T':[ nTA/totalrefT , nTC/totalrefT , nTG/totalrefT ,             0 ] }

ffactor = map(lambda x: 1.0/sum(x),GTR_nonCpG.values())
for ki,key in enumerate(GTR_nonCpG.keys()):
  for i in range(4):
    GTR_nonCpG[key][i] = GTR_nonCpG[key][i] * ffactor[ki]

##SUBSTITUTION MATRIX CpG
mtotalrefC = (totalCpG-(mCA+mCG+mCT+mGA+mGC+mGT))/2.0+mCA+mCG+mCT
mtotalrefG = (totalCpG-(mCA+mCG+mCT+mGA+mGC+mGT))/2.0+mGA+mGC+mGT

GTR_CpG = {   'C':[ mCA/mtotalrefC ,              0 , mCG/mtotalrefC , mCT/mtotalrefC ],
              'G':[ mGA/mtotalrefG , mGC/mtotalrefG ,              0 , mGT/mtotalrefG ] }

ffactor = map(lambda x: 1.0/sum(x),GTR_CpG.values())
for ki,key in enumerate(GTR_CpG.keys()):
  for i in range(4):
    GTR_CpG[key][i] = GTR_CpG[key][i] * ffactor[ki]

# INITIATE INSERT NORMALIZED LIKELIHOODS
inserts = []
tinserts = float(sum(insertionsizes.values()))
sys.stderr.write("TOTAL INSERTS FROM EPO:\t%d\n"%tinserts)
to_sort = insertionsizes.keys()
to_sort.sort()
sys.stderr.write("Length\tCount\tnFreq\n")
for key in to_sort:
  value = insertionsizes[key]/tinserts
  inserts.append((key,value))
del insertionsizes
ginserts = tinserts/sum_obs
sys.stderr.write("INSERTIONS: gfreq=%.8f, %s...\n"%(ginserts,str(inserts)[:120]))

# INITIATE DELETION NORMALIZED LIKELIHOODS
deletions = []
tdeletions = float(sum(deletionsizes.values()))
sys.stderr.write("TOTAL DELETIONS FROM EPO:\t%d\n"%tdeletions)
to_sort = deletionsizes.keys()
to_sort.sort()
sys.stderr.write("Length\tCount\tnFreq\n")
for key in to_sort:
  value = deletionsizes[key]/tdeletions
  #if key <= 10: sys.stderr.write("%d\t%d\t%.8f\n"%(key, deletionsizes[key], value))
  deletions.append((key,value))
del deletionsizes
gdeletions = tdeletions/sum_obs
sys.stderr.write("DELETIONS: gfreq=%.8f, %s...\n"%(gdeletions,str(deletions)[:120]))

#OPEN FASTA COPY WITH MMAP
sys.stderr.write('Doing mmap of genome file...\n')
f = open(options.infile, "r")
cmap = mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ)
cmap.flush()

sys.stderr.write('Iterating over genome...\n')

VCFheader = """##fileformat=VCFv4.1
##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context.">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\torig\tmutant"""
output.write(VCFheader+"\n")

chrom,length,sblock,bline,cline = fastaindex

if total > 1000: y = mut/float(total)
else: y=gmut
if totalCpG > 1000: yCpG = mutCpG/float(totalCpG)
else: yCpG=gmutCpG

#dmut, gmut, gdeletions (deletions), ginserts (inserts)
lmut = y/gmut * dmut
lmutGpG = yCpG/gmutCpG * dmutCpG
ldeletion = y/gmut * gdeletions/gmut * dmutindel
linsertion = y/gmut * ginserts/gmut * dmutindel

lbase,base,nbase = None,None,None
for pos in xrange(2,length-1):
  if lbase != None:
    lbase,base,nbase = base,nbase,None
    npos = (pos // bline)*cline+(pos % bline)
    cmap.seek(sblock+npos)
    nbase = cmap.read_byte().upper()
  else:
    npos = ((pos-2) // bline)*cline+((pos-2) % bline)
    cmap.seek(sblock+npos)
    lbase = cmap.read_byte().upper()
    npos = ((pos-1) // bline)*cline+((pos-1) % bline)
    cmap.seek(sblock+npos)
    base = cmap.read_byte().upper()
    npos = (pos // bline)*cline+(pos % bline)
    cmap.seek(sblock+npos)
    nbase = cmap.read_byte().upper()

  # GTR-CpG mutations
  ismut = random.random()
  if (base == "C" and nbase == "G") or (base == "G" and lbase == "C"):
    if ismut <= lmutGpG:
      whichmut = random.random()
      for i,alt in enumerate('ACGT'):
        if alt != base:
          whichmut-=GTR_CpG[base][i]
          if whichmut < 0:
            output.write("%s\t%d\t.\t%s\t%s\t.\t.\tCpG\tGT\t0\t1\n"%(chrom,pos,base,alt))
            break
  # GTR mutations
  elif base in "ACGT":
    if ismut <= lmut:
      whichmut = random.random()
      for i,alt in enumerate('ACGT'):
        if alt != base:
          whichmut-=GTR_nonCpG[base][i]
          if whichmut < 0:
            output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t0\t1\n"%(chrom,pos,base,alt))
            break
  else:
    #sys.stderr.write("Base %s on %s %d can not be mutated. Skipping.\n"%(base,chrom,pos))
    continue

  # Indels
  isindel = random.random()
  if isindel <= ldeletion:
    npos = ((pos-1) // bline)*cline+((pos-1) % bline)
    cmap.seek(sblock+npos)
    base = cmap.read_byte().upper()
    if base == "N": continue
    indellength = random.random()
    for dellength,value in deletions:
      indellength-=value
      if indellength <= 0:
        npos = (pos // bline)*cline+(pos % bline)
        cmap.seek(sblock+npos)
        delseq = ''
        for i in xrange(dellength):
          if pos+i < length:
            npos = ((pos+i) // bline)*cline+((pos+i) % bline)
            cmap.seek(sblock+npos)
            delseq += cmap.read_byte().upper()
          else:
            break
        if len(delseq) == dellength:
          output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t0\t1\n"%(chrom,pos,base+delseq,base))
        break

  isindel = random.random()
  if isindel <= linsertion:
    npos = ((pos-1) // bline)*cline+((pos-1) % bline)
    cmap.seek(sblock+npos)
    base = cmap.read_byte().upper()
    if base == "N": continue
    indellength = random.random()
    for inslength,value in inserts:
      indellength-=value
      if indellength <= 0:
        insseq = ""
        for s in range(inslength):
          whichbase = random.random()
          for i,newbase in enumerate('ACGT'):
            whichbase -= bfreqs[i]
            if whichbase <= 0:
              insseq += newbase
              break
        output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t0\t1\n"%(chrom,pos,base,base+insseq))
        break

cmap.close()
f.close()
output.close()

# Write a fasta with the mutations using GATK
if options.fasta_out == True:
  sys.stderr.write("Creating fasta with mutations\n")
  try:
    gatk_command = java_location + " -Xmx300M -jar " + gatk_location + " -R " + options.infile + " -T FastaAlternateReferenceMaker -o " + options.outfile + ".fa --variant " + options.outfile + ".vcf"
    retcode = subprocess.call(gatk_command, shell=True)
    if retcode < 0:
        print >>sys.stderr, "GATK was terminated by signal", -retcode
  except OSError as e:
    print >>sys.stderr, "Execution failed:", e
