
library(SpatialSimEx) # This is a special R-package created for this example. 
library(mvtnorm)
library(fields)
library(Matrix)
library(matrixcalc)
library(nlme)
library(corpcor)
library(acebayes)
library(doParallel)

source("combine_ute.R")
source("Exp_crit.R")
source("Laplace_approx.R")
source("log_posterior.R")
source("trace_mat.R")
source("Predict_uteApp.R")

cl <- makeCluster(25,outfile="cl_out.txt")
registerDoParallel(cl)

N <- 5000 # number of particles
Scenario <- 1 # {1,2,3}
if(Scenario==1){
  Sp.range <- 0.2
}else if(Scenario==2){
  Sp.range <- 0.5
}else{
  Sp.range <- 0.8
}

mu_prior1 <- c(5,-2.8,8,log(1.2),log(.7/Sp.range),log(Sp.range))  # prior means for model 1
Sigma_prior1 <- c(4,4,4,0.25,0.25,0.25)   # prior variances for model 1

mu_prior2 <- c(3.8,-1.5,-1.2,log(.6/Sp.range),log(Sp.range))   # prior means for model 2
Sigma_prior2 <- c(0.125,0.125,0.125,0.125,0.125)   # prior variances for model 2

mu_Ltau <- log(.7/.3) # prior mean for logit(tau)
Sigma_Ltau <- 0.25    # prior variance for logit(tau)

mu_prior <- c(mu_prior1,mu_prior2,mu_Ltau)
Sigma_prior <- diag(c(Sigma_prior1,Sigma_prior2,Sigma_Ltau))

prior_b10 <- rnorm(N,mu_prior[1],sqrt(Sigma_prior[1,1]))
prior_b11 <- rnorm(N,mu_prior[2],sqrt(Sigma_prior[2,2]))
prior_b12 <- rnorm(N,mu_prior[3],sqrt(Sigma_prior[3,3]))
prior_varY <- rnorm(N,mu_prior[4],sqrt(Sigma_prior[4,4]))
prior_r11.r12 <- rnorm(N,mu_prior[5],sqrt(Sigma_prior[5,5]))
prior_r12 <- rnorm(N,mu_prior[6],sqrt(Sigma_prior[6,6]))

prior_b20 <- rnorm(N,mu_prior[7],sqrt(Sigma_prior[7,7]))
prior_b21 <- rnorm(N,mu_prior[8],sqrt(Sigma_prior[8,8]))
prior_b22 <- rnorm(N,mu_prior[9],sqrt(Sigma_prior[9,9]))
prior_r21.r22 <- rnorm(N,mu_prior[10],sqrt(Sigma_prior[10,10]))
prior_r22 <- rnorm(N,mu_prior[11],sqrt(Sigma_prior[11,11]))

logit_tau <- rnorm(N,mu_prior[12],sqrt(Sigma_prior[12,12]))

# Initial particle set
prior <- data.frame(prior_b10,prior_b11,prior_b12,prior_varY,prior_r11.r12,prior_r12,prior_b20,prior_b21,prior_b22,prior_r21.r22,prior_r22,logit_tau)

iSigma_prior <- solve(Sigma_prior)
log_det_prior <- log(det(Sigma_prior))

num_par <- ncol(prior)  # total number of parameters in the Copula model
r01 <- 0.001            # nugget effect parameter for model 1
r02 <- 0.001            # nugget effect parameter for model 2

# Prediction locations
Unsamp_Lx <- rep(seq(0,1,.25),5)
Unsamp_Ly <- rep(seq(0,1,.25),each=5)

X1 <- ((Unsamp_Ly+.1)*5)-((Unsamp_Lx+.2)*5) # values of the 1st covariate at the prediction locations
X2 <- (Unsamp_Ly+.1)*(Unsamp_Lx+.2)         # values of the 2nd covariate at the prediction locations

site_no <- 1:25 # site numbers

Unsamp_X <- data.frame(site_no,Unsamp_Lx,Unsamp_Ly,X1,X2)
dist_mat <- rdist(data.frame(Unsamp_Lx,Unsamp_Ly))

################ Optimal design selection using ACE algorithm ###############

d_no <- 5 # number of design points
dint <- matrix(runif(d_no*2),ncol=2) # initial design

set.seed(101)
tot_samples <- 1000

B1 <- 2000
theta_int <- rmvnorm(n=B1,mean=mu_prior, sigma = Sigma_prior)

B2 <- 1000
Z1 <- matrix(rnorm(tot_samples*B2), tot_samples)

Unif_set <- matrix(runif(tot_samples*2),ncol=2)
Z_res <- matrix(rnorm(tot_samples*2),ncol=2)

opt_data <- ace(exp.crit, start.d=dint, B=c(2000,500), Q = 20, N1 = 10, N2 = 0, lower=0, upper = 1)

save(opt_data, file = "out.RData")



