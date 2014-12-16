#
# Libraries
#
library(coda)
library(rjags)

# Model based on elements of code for exercise 9.2 in 
# 'Doing Bayesian Data Analysis', 1st edition, Kruschke 2011

#
# Constants
#
data_location <- "~/Documents/PhD/hsd_locus/mapping/"

five_prime_input <- paste(data_location, "jags_5prime_input.txt", sep="")
three_prime_input <- paste(data_location, "ivr_mapped_alleles.txt", sep="")

# Possible outcomes
alleles = c("A", "B", "C", "D", "E", "F")

#
# JAGS chain parameters
#
parameters = c("mu", "kappa", "theta", "a", "b") # Parameters to output posterior distributions

# MCMC params
adapt_steps = 500 
burn_in_steps = 2000
num_chains = 3
num_save_steps = 50000
thin_steps = 1

# Steps per chain
num_iterations = ceiling((num_save_steps * thin_steps ) / num_chains)

#
# Functions
#
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
# Run model for first allele
#

#
# Hierarchical model spec
#
jags_model1_spec = "
# JAGS model specification
model {
  # For each sample, likelihood and prior
  for (sample_index in 1:num_samples)
  {
    # likelihood
    # Number of reads that map to one allele is binomial (i.e. 1.1 = success, 1.2 = failure)
    y[sample_index] ~ dbin(theta[sample_index], N[sample_index]) 
    
    # Beta prior for proportion of population with allele in each sample
    theta[sample_index] ~ dbeta(a[tissue[sample_index]], b[tissue[sample_index]])T(0.0001,0.9999)
  }

  # For each tissue (blood or csf) hyperpriors
  for (tissue_index in 1:num_tissues)
  {
    # Convert a and b in beta prior, to mu and kappa 
    # mu = mean allele for tissue, 
    # kappa = how closely does sequence represent tissue - constant across all tissue types
    a[tissue_index] <- mu[tissue_index] * kappa
    b[tissue_index] <- (1 - mu[tissue_index]) * kappa

    # hyperpriors for mu (beta dist) 
    mu[tissue_index] ~ dbeta(Amu, Bmu)
  }

  # Kappa hyperprior (gamma dist - shape and rate)
  kappa ~ dgamma(Skappa, Rkappa)

  # Top level constants for hyperpriors
  # Beta dist for mu - estimated from Manso et al 2014 fig 4h
  # 84 mice, 5% representation of allele 1.2 in blood
  Amu <- 80
  Bmu <- 4

  # Gamma dist for kappa. First convert mean and sd to shape and rate
  Skappa <- pow(meanGamma,2)/pow(sdGamma,2)
  Rkappa <- meanGamma/pow(sdGamma,2)

  meanGamma <- 10
  sdGamma <- 10
}
"
writeLines(jags_model1_spec,con="model1.txt")
#
# Read in data
#
five_prime_reads <- read.delim(five_prime_input)

# Convert to a list for use with JAGS
five_prime_data = list(num_tissues = length(unique(five_prime_reads$Tissue)), num_samples = length(five_prime_reads$TotalReads), tissue = five_prime_reads$Tissue, N = five_prime_reads$TotalReads, y = five_prime_reads$AReads)

#
# Run the model
#

# Create and adapt
cat("Running first model\n\n")
jags_model1 = jags.model("model1.txt", data=five_prime_data, n.chains=num_chains, n.adapt=adapt_steps)

# Burn-in
cat("MCMC burn in iterations...\n")
update(jags_model1, n.iter=burn_in_steps)

# Converged chain
cat( "Sampling iterations...\n" )
coda_samples1 = coda.samples(jags_model1, variable.names=parameters, n.iter=num_iterations, thin=thin_steps)

# Save the chain
save(coda_samples1, file="chain1.Rdata")

#
# Use first model posteriors to produce data for second model
#
cat("Converting data\n\n")
mcmc_chain = as.matrix(coda_samples)

# Output from ivr-typer
# Columns: sample, 1.1 reads, 1.2 reads, 2.1 reads, 2.2 reads, 2.3 reads, 2.1 reads, 2.2 reads, 2.3 reads
# Where the last 6 fields are 3 mapping to 1.1, and 3 mapping to 1.2 respectively
#
# Need to remove samples with no reads mapping in the first allele before processing
# awk '$2+$3!=0 {print $0}' ivr_mapped_alleles.txt > jags_3prime_input.txt
hsd_mapping <- read.delim(three_prime_input, header=F)

# Copy five prime data structure
three_prime_reads <- five_prime_reads[c("TotalReads", "Tissue")]
N <- NULL

# For each sample in the first chain
for (i in 1:nrow(three_prime_reads)) {
  # Extract the column of the matrix
  col_name <- paste("theta[",i,"]",sep="")
  theta = mcmc_chain[,col_name]
  
  # Sums of reads mapping to 1.1 and 1.2 respectively
  # This gives a distribution for sample size.
  reads <- hsd_mapping[i,seq(from=4, to=9)]
  Ndist <- theta*(sum(reads[c("V4", "V5", "V6")])) + (1-theta)*(sum(reads[c("V7", "V8", "V9")]))
  # Take a sample from it
  N[i] <- round(sample(Ndist,1))
  
  # Number of reads mapping to 1.1, sampling from posterior for theta
  allele1 <- sum(sample(theta, N[i], replace=TRUE) > 0.5)
  
  # Set up weighted sample for counts of alleles.
  weights1 = vet_weights(reads[c("V4","V5","V6")])
  weights2 = vet_weights(reads[c("V7","V8","V9")])
  # Output is a table of counts for the six possible alleles
  sampled_alleles <- table(c(sample(c("A","B","E"),allele1,replace=TRUE,prob=weights1), sample(c("D","C","F"),N[i]-allele1,replace=TRUE,prob=weights2)))
  
  # Convert this table into a data frame (badly)
  for (j in 1:length(alleles))
  {
    if (is.na(sampled_alleles[alleles[j]]))
    {
      three_prime_reads[i,alleles[j]] = 0
    }
    else
    {
      three_prime_reads[i,alleles[j]] = sampled_alleles[alleles[j]]
    }
  }
  # Free memory where possible
  rm(theta)
}

rm(mcmc_chain)

three_prime_reads$TotalReads <- N

# Save the converted data
save(three_prime_reads, file="three_prime_reads.Rdata")

#
# Run second model
#
jags_model2_spec = "
# JAGS model specification
model {
  # For each sample, likelihood and prior
  for (sample_index in 1:num_samples)
  {
    # likelihood
    # Number of reads that map to each allele is multinomial. 
    # y and pi are vectors of length 6
    y[sample_index] ~ dmulti(pi[sample_index], N[sample_index]) 
    
    # Beta prior for proportion of population with allele in each sample
    pi[sample_index] ~ ddirch(alpha[tissue[sample_index]])T(0.0001,0.9999)
  }

  # For each tissue (blood or csf) hyperpriors
  for (tissue_index in 1:num_tissues)
  {
    # Convert a and b in beta prior, to mu and kappa 
    # mu = mean allele for tissue, 
    # kappa = how closely does sequence represent tissue - constant across all tissue types
    alpha[tissue_index] <- mu[tissue_index] * kappa

    # hyperpriors for mu (beta dist) 
    mu[tissue_index] ~ ddirch(AlphaMu)
  }

  # Kappa hyperprior (gamma dist - shape and rate)
  kappa ~ dgamma(Skappa, Rkappa)

  # Top level constants for hyperpriors
  # Beta dist for mu - estimated from Manso et al 2014 fig 4h
  # 84 mice. A: 20; B: 10; C: 1; D: 3; E: 65; F: 1
  AlphaMu <- c(20, 10, 1, 3, 65, 1)

  # Gamma dist for kappa. First convert mean and sd to shape and rate
  Skappa <- pow(meanGamma,2)/pow(sdGamma,2)
  Rkappa <- meanGamma/pow(sdGamma,2)

  meanGamma <- 10
  sdGamma <- 10
}
"
writeLines(jags_model2_spec,con="model2.txt")

# Convert to a list for use with JAGS
three_prime_data = list(num_tissues = length(unique(three_prime_reads$Tissue)), num_samples = length(three_prime_reads$TotalReads), tissue = three_prime_reads$Tissue, N = three_prime_reads$TotalReads, y = three_prime_reads[,alleles])

#
# JAGS chain parameters
#
parameters = c("mu", "kappa", "pi", "alpha") # Parameters to output posterior distributions

#
# Run the model
#

# Create and adapt
cat("Running second model\n\n")
jags_model2 = jags.model("model2.txt", data=three_prime_data, n.chains=num_chains, n.adapt=adapt_steps)

# Burn-in
cat("MCMC burn in iterations...\n")
update(jags_model2, n.iter=burn_in_steps)

# Converged chain
cat( "Sampling iterations...\n" )
coda_samples2 = coda.samples(jags_model2, variable.names=parameters, n.iter=num_iterations, thin=thin_steps)

# Save the chain
save(coda_samples2, file="chain2.Rdata")
