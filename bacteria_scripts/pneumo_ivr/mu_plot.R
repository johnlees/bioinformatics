# Plot output from mcmc chain, for my per sample

# Start by converting the relevant parameters for mu
mu_diff = as.matrix(coda_samples[,c(6,7),drop=FALSE])

# Significant difference between tissues?
hist(mu_diff[,1]-mu_diff[,2])

# Plot a histogram of both together
tissue_header <- c(rep("csf",length(df$csf)), rep("blood",length(df$blood)))

mu_df<-as.data.frame(c(mu_diff[,1], mu_diff[,2]))
mu_df$tissue <- tissue_header
names(mu_df) <- c("mu", "tissue")

ggplot(mu_df) + geom_density(aes(x=mu, fill=tissue),alpha=0.5)