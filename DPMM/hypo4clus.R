library("MASS")
# simulate 4 clusters of data
set.seed(11)
n <- 60
m1 <- c(1.5, 1.5)   # upper right x, y means
S1 <- matrix(c(0.3, 0.05, 0.05, 0.3), ncol = 2)
clus1 <- mvrnorm(n = n, mu = m1, Sigma = S1)
m2 <- c(1.5, -1.5)  # lower right
S2 <- matrix(c(0.5, -0.08, -0.08, 0.2), ncol = 2)
clus2 <- mvrnorm(n = n, mu = m2, Sigma = S2)
m3 <- c(-1.5, 1.5)  # upper left
S3 <- matrix(c(0.1, 0.03, 0.03, 0.1), ncol = 2)
clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)
m4 <- c(-1.5, -1.5) # lower left
S4 <- matrix(c(0.8, 0.50, 0.50, 0.8), ncol = 2)
clus4 <- mvrnorm(n = n, mu = m4, Sigma = S4)
datc <- rbind(clus1, clus2, clus3, clus4)

# import the crp_gibbs() function 

# 
# Dirichlet Process Mixture Model example. 
# Last modified on Apr 19, 2018 by Yuelin Li from an example given by 
# Tamara Broderick, whose original source code is available at: 
# https://raw.githubusercontent.com/tbroderick/bnp_tutorial
#      /2016mlss/ex7_sampler.R 
###### Below are headers in the original source code file ###### 
# A Gibbs sampler for a CRP Gaussian mixture model 
# Algorithm 3 in Neal 2000 
# Copyright (C) 2015, Tamara Broderick 
# www.tamarabroderick.com 
# 
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <-ttp://www.gnu.org/licenses/>. 
# useful 	for working with multivariate normal distributions 
#
##### Below are modified by YL. Variable names are changed from 
# the original to match the best we can the notation in the tutorial. 
crp_gibbs <- function(data,alpha=0.01,mu0,sigma0,sigma_y,c_init,maxIters=1000)
{ 
# data: an N by D matrix of data points 
# A small alpha encourages fewer clusters. 
#    alpha <- 0.01 
# sigma_y: measurement error of data y, assumed known, same for all clusters.
#    sigma_y <- diag(data_dim) * 1
# mu0, sigma0: prior mean and variance around unknown mean mu
#    mu0 <- matrix(rep(0, data_dim), ncol = data_dim, byrow = TRUE) 
#    sigma0 <- diag(data_dim) * 3^2
# c_init: initial assignments of data points to clusters
# 
require(mvtnorm) 
# dimension of the data points 
data_dim <- ncol(data) 
# number of data points 
N <- nrow(data) 
##### 
# Priors 
##### 
# prior precision on the unknown mean mu.
tau0 <- solve(sigma0) # prior precision on $mu$, inverse of prior covariance 
# cluster-specific precision, assumed known, all clusters are assumed to 
# share identical measurement error of y ~ N(mu, sigma_y).
tau_y <- solve(sigma_y) 
# initialize the CRP Gibbs sampler 
z <- c_init                # initial cluster membership assignments  
n_k <- as.vector(table(z)) # initial data counts at each cluster, Eq (4). 
Nclust <- length(n_k)      # initial number of clusters 
##
# Chinese Restaurant Process (CRP) Gibbs sampler begins 
##
res <- matrix(NA, nrow = N, ncol = maxIters) # cluster membership storage 
pb <- txtProgressBar(min = 0, max = maxIters, style = 3) 
for(iter in 1:maxIters) { # maxIters also prevents endless loops 
for(n in 1:N) {    # one data point (customer) at a time 
# when nth customer enters the Chinese restaurant, we need to first 
# un-assign his/her initial cluster membership, then use Eq (12) 
# [already occupied table] and (13) [new table] to calculate the 
# updated probability of table assignment. 
c_i <- z[n] # what is the nth person's table assignment? 
n_k[c_i] <- n_k[c_i] - 1 # remove the nth person from table 
# if the table becomes empty when the nth person is removed, 
# then that table/cluster is removed. 
if( n_k[c_i] == 0 ) 
{ 
  n_k[c_i] <- n_k[Nclust]   # last cluster to replace this empty cluster 
  loc_z <- ( z == Nclust )  # who are in the last cluster? 
  z[loc_z] <- c_i           # move them up to fill just emptied cluster 
  n_k <- n_k[ -Nclust ]     # take out the last cluster, now empty 
  Nclust <- Nclust - 1      # decrease total number of clusters by 1 
} 
z[n] <- -1 # ensures z[n] doesnÕt get counted as a cluster ##### 
##
# Now we are ready to update table assignment by Eqs (12) and (13). 
##
# log probabilities for the clusters, add previously unoccupied table 
logp <- rep( NA, Nclust + 1 ) 
# loop over already occupied tables 1:J and calculate pr as per Eq (13). 
for( c_i in 1:Nclust ) { 
  tau_p <- tau0 + n_k[c_i] * tau_y  # cluster precision as per Eq (4)
  sig_p <- solve(tau_p)             # cluster variance, inverse of tau_c
  # find all of the points in this cluster 
  loc_z <- which(z == c_i) 
  # sum all the points in this cluster 
  if(length(loc_z) > 1) { 
    sum_data <-  colSums(data[z == c_i, ]) } 
  else { 
    sum_data <-  data[z == c_i, ] 
  } 
# 
# We need to use the predictive distribution of each already 
# occupied table to predict the next customer sitting there. 
# 
# Use Eq (4) to find the conditional posterior distribution for 
# the cluster means, then use the predictive distribution of 
# y_j in Eq (12) to predict new data value c_i from c-i. 
# 
mean_p <-  sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0)) 
logp[c_i] <- log(n_k[c_i]) + 
   dmvnorm(data[n,], mean = mean_p, sigma = sig_p + sigma_y, log = TRUE) } 
# 
# We are done looping over already occupied tables. Next, we use 
# Eq (13) to calcualte the log probability of a previously 
# unoccupied, "new" table. Essentially, it is the prior predicitive 
# distribution of the DP.
# 
logp[ Nclust+1 ] <- log(alpha) + 
    dmvnorm(data[n,], mean = mu0, sigma = sigma0 + sigma_y, log = TRUE) 
# transform unnormalized log probabilities into probabilities 
max_logp <- max(logp) 
logp <- logp - max_logp 
loc_probs <- exp(logp) 
loc_probs <- loc_probs / sum(loc_probs) 
# draw a sample of which cluster this new customer should belong to 
newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs) 
# spawn a new cluster if necessary 
if(newz == Nclust + 1) { 
  n_k <- c(n_k, 0) 
  Nclust <- Nclust + 1 
  } 
z[n] <- newz 
n_k[newz] <- n_k[newz] + 1 # update the cluster n_k 
} 
setTxtProgressBar(pb, iter) # update text progress bar after each iter 
res[, iter] <- z     # cluster membership of N observations
} 
close(pb)            # close text progress bar 
invisible(res)       # return results, N by maxIters matrix 
}
#####  end of crp_gibbs() #####

# run crp_gibbs() 
alpha <- 0.01
mu0 <- matrix(rep(0, 2), ncol = 2, byrow = TRUE)
sigma0 <- diag(2) * 3^2
sigma_y <- diag(2) * 1
c_init <- rep(1, nrow(datc))
results <- crp_gibbs(data = datc, alpha = alpha,
 mu0 = mu0, sigma0 = sigma0, sigma_y = sigma_y, 
      c_init = rep(1, nrow(datc)))
##### > # Gibbs sampling iterations are saved in
# ‘results’. Each row contains the latent
# class membership distribution of one
# person over MCMC iterations. A person
# is deemed a member of c=2 if, e.g., s/he
# gets assigned most often to c2. So we
# tabulate the frequency counts and
# whichever cluster with highest frequency
# is that person’s latent class.
#####
tab <- apply( results, 1, FUN = function(x) { tab <- table(x)
    ans <- names(tab[which.max(tab)])
    return(ans) } )

print(table(tab))





