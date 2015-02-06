#
# A Bayesian hierarchical model for hsdS allele in S. Pneumo, based on the D39 reference
# Estimates the allele distribution both within a single sample, and in the tissue type
# it was taken from (e.g. blood, csf, nasopharynx)
#
# Author: John Lees
# 2014-2015
#
# Model based on elements of code for exercise 9.2 in
# 'Doing Bayesian Data Analysis', 1st edition, Kruschke 2011
#

#
# Libraries
#
library(rstan)
library(parallel)

#
# Constants
#

# Data input
data_location <- "~/Documents/PhD/hsd_locus/mapping"

five_prime_input <- paste(data_location, "jags_5prime_input.txt", sep="/")
three_prime_input <- paste(data_location, "jags_3prime_input.txt", sep="/")

# Possible outcomes
alleles = c("A", "B", "C", "D", "E", "F")

# Seed for rng - set to make reproducible
# Using rjags alone paralellises with plyr automatically, but seed is the same and
# all chains are therefore the same.
# Instead parallelise with snow and ClusterMap
rng_seed = 1
set.seed(rng_seed)

# default MCMC params
num_chains = 3

#
# Functions
#

# Assume column vectors that sum to zero are equally weighted
vet_weights <- function(x) {
  if (sum(x) == 0)
  {
    re_weighted <- rep("1",times=3)
  }
  else
  {
    re_weighted <- x
  }
  return(re_weighted)
}

#
# Main
#

#
# Run model for first allele
#

# Hierarchical model spec
stan_model1_spec = "
# STAN model specification
data {
  int<lower=0> num_tissues;
  int<lower=0> num_samples;
  int<lower=0> tissue[num_samples];
  int<lower=0> N[num_samples];
  int<lower=0> y[num_samples];
}
parameters {
  real<lower=0.0001,upper=0.9999> theta[num_samples];
  real<lower=0> kappa;
  real<lower=0.0001,upper=0.9999> mu[num_tissues];
}
transformed parameters {
  real<lower=0> a[num_tissues];
  real<lower=0> b[num_tissues];
  for (j in 1:num_tissues) {
    a[j] <- mu[j] * kappa;
    b[j] <- (1-mu[j]) * kappa;
  }
}
model {
  kappa ~ gamma(1, 0.1);

  for (j in 1:num_tissues) {
    mu[j] ~ beta(80, 4);
  }

  for (j in 1:num_samples) {
    theta[j] ~ beta(a[tissue[j]], b[tissue[j]]);
    y[j] ~ binomial(N[j], theta[j]);
  }
}
"
writeLines(stan_model1_spec,con="stan_model1.txt")
#
# Read in data
#
five_prime_reads <- read.delim(five_prime_input)

# Convert to a list for use with JAGS
five_prime_data = list(num_tissues = length(unique(five_prime_reads$Tissue)),
  num_samples = length(five_prime_reads$TotalReads), tissue = five_prime_reads$Tissue,
  N = five_prime_reads$TotalReads, y = five_prime_reads$AReads)

#
# Run the model
#
fit1 <- stan(file = 'stan_model1.txt', data = five_prime_data, 
            iter = 1000, chains = 0)
sflist <- 
    mclapply(1:4, mc.cores = 2, 
             function(i) stan(fit = fit1, data = five_prime_data, 
                              seed = rng_seed, 
                              chains = 1, chain_id = i, 
                              refresh = -1))

fit <- sflist2stanfit(sflist)

pdf("model1.pdf")
plot(fit)
dev.off()

pdf("model1_trace.pdf")
traceplot(fit)
dev.off()

capture.output(print(fit,max=20000),file="param_summary.txt")
