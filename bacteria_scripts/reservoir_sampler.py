#!/usr/bin/env python
# -*- coding: ASCII -*-

#
# Based on http://en.wikipedia.org/wiki/Reservoir_sampling
#

import random
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n","--samples", dest="samples", help="Number of samples to subsample",type="int",default=100)
parser.add_option("-f","--file", dest="file", help="File to sample from")
parser.add_option("-s","--seed", dest="seed", help="Seed for random number generator",type="int",default=1)
(options, args) = parser.parse_args()
 
# Force the value of the seed so the results are repeatable
random.seed(options.seed)
 
infile = open(options.file)
samples = []

for index, line in enumerate(infile):
        # Generate the reservoir
        if index < options.samples:
                samples.append(line.rstrip('\n'))
        else:                  
                # Randomly replace elements in the reservoir
                # with a decreasing probability.             
                # Choose an integer between 0 and index (inclusive)               
                r = random.randint(0, index)               
                if r < options.samples:                       
                        samples[r] = line.rstrip('\n')

infile.close()

print "\n".join(samples)