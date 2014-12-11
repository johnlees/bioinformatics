library(rjags)

# Based on code for exercise 9.2 in 'Doing Bayesian Data Analysis',
# 1st edition, Kruschke 2011

#
# Hierarchical model spec
#
jagsModelSpec = "
# JAGS model specification
model {
  # For each sample, likelihood and prior
  for (sample_index in 1:num_samples)
  {
    # likelihood
    # Number of reads that map to one allele is binomial (i.e. 1.1 = success, 1.2 = failure)
    y[sample_index] ~ dbin(theta[sample_index], N[sample_index]) 
    
    # Beta prior for proportion of population with allele in each sample
    theta[sample_index] ~ dbeta(a[tissue[sample_index]], b[tissue[sample_index]])I(0.001,0.999)
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
writeLines(modelString,con="model.txt")
#
# Read in data
#
5primeReads <- read.table("5prime_jags.txt", header=T)

# Convert to a list for use with JAGS
5primeData = list(num_tissues = length(unique(5primeReads$Tissue)), num_samples = length(5primeReads$TotalReads), tissue = 5primeReads$Tissue, N = 5primeReads$TotalReads, y = 5primeReads$AReads)

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
# Run the model
#

# Create and adapt
jags_model = jags.model("model.txt" , data=5primeData, n.chains=num_chains, n.adapt=adapt_steps)

# Burn-in
cat("MCMC burn in iterations...\n")
update(jags_model , n.iter=burn_in_steps)

# Converged chain
cat( "Sampling iterations...\n" )
coda_samples = coda.samples(jags_model, variable.names=parameters, n.iter=nIter, thin=thin_steps)

# resulting codaSamples object has these indices:
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

# TODO: should check chain convergence and autocorrelation

# Convert to matrix, though note each chain is concatenated
mcmc_chains = as.matrix( coda_samples )

# Save the results
mu = NULL
for (i in 1:num_tissues)
{
  mu = rbind(mu, mcmc_chains(paste(",mu[",i,"]",sep="")))
}

kappa = mcmc_chains[,"kappa"]
save(mu, kappa, file="jags_result.Rdata")

