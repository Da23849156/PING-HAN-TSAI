# DPpackage is no longer available on CRAN because it is not actively 
# maintained.  However, it can be found on GitHub.  It needs compiler 
# tools to run.  On Windows 10, you need to get Rtools from 
# https://cran.r-project.org/bin/windows/Rtools/.   Once Rtools are
# installed on Windows 10, you can download and then compile
# DPpackage using this command:
# install.packages("C:\\DPpackage_1.1-7.tar.gz", repos=NULL,type="source")
#
# On Mac, Xcode is needed to compile DPpackage source code.
# (see https://clanfear.github.io/CSSS508/docs/compiling.html).
#
# On my old Macbook Pro running High Sierra (OS 10.13.6), I had to 
# manually download and install 3 tools for DPpackage to compile.
# 1. Xcode 10.1 (an older version) from Apple Developer site
#    https://developer.apple.com/download/more/
#    (note: newer versions of Xcode produce an error after unzip)
# 2. "Command Line Tools (macOS 10.13) for Xcode 10.1", and 
# 3. gfortran for High Sierra from a different website:
#    https://github.com/fxcoudert/gfortran-for-macOS/releases
#

insted_pkgs <- installed.packages()
# install needed packages if they are not found
if ( (! 'devtools' %in% insted_pkgs) ) install.packages("devtools")
# A 3rd way to install DPpackage is to get it from GitHub:
if ( (! 'DPpackage' %in% insted_pkgs) ) 
	devtools::install_github("konkam/DPpackage")

# This script is adapted from the example in the DPpackage documentation.
# https://users.soe.ucsc.edu/~draper/DPpackage.pdf
library("DPpackage")
data(galaxy)
prior2 <- list(alpha=1,          # controls number of clusters in DP prior
    m1=rep(0,1),                 # G0 mean (spawning new mixture component)
    psiinv2=solve(diag(0.5,1)),  # hyperprior on precision of G0
    nu1=4,nu2=4,tau1=1,tau2=100) # nu1: mean of precision of G0
                     # nu2: hyperprior on mean of psi1
		     # tau1, tau2: gamma prior on k0, scale of Normal part
		     #             of G0.

data(galaxy)
galaxy <- data.frame(galaxy, speeds = galaxy$speed/1000)
attach(galaxy)
state <- NULL        # Initial state
# MCMC parameters
nburn <- 1000
nsave <- 10000
nskip <- 10
ndisplay <- 100
mcmc <- list(nburn = nburn, nsave = nsave, nskip = nskip, ndisplay = ndisplay)
# Example of Prior information 2
# Fixing alpha and m1
prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),nu1=4,nu2=4,
	       tau1=1,tau2=100)
# Example of Prior information 3, Fixing only alpha
prior3 <- list(alpha=1,m2=rep(0,1),s2=diag(100000,1),
	       psiinv2=solve(diag(0.5,1)),
	       nu1=4,nu2=4,tau1=1,tau2=100)
# Example of Prior information 4, Everything is random
prior4 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
	       psiinv2=solve(diag(0.5,1)),
	       nu1=4,nu2=4,tau1=1,tau2=100)
# Fit the models
fit1.2 <- DPdensity(y=speeds,prior=prior2,mcmc=mcmc,state=state,status=TRUE)
fit1.3 <- DPdensity(y=speeds,prior=prior3,mcmc=mcmc,state=state,status=TRUE)
fit1.4 <- DPdensity(y=speeds,prior=prior4,mcmc=mcmc,state=state,status=TRUE)
# Posterior means
fit1.2
fit1.3
fit1.4
# Plot distributions of model parameters.  Note density of categories.
# Thus, fit1.2$dens is a weighted density over all BNP categories.
# I haven't been able to find each observation's latent class membership.
plot(fit1.2, output="param")
# Plot the estimated density
hist(galaxy$speeds, probability = T, breaks = seq(5, 35, by =1), col = "grey80", border = "grey60")
lines(cbind(fit1.2$x1,fit1.2$dens), lwd = 3)

