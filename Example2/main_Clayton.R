
library(AirQualityRcpp) #  This is a special R-package created for this example.
library(mvtnorm)
library(fields)
library(Matrix)
library(matrixcalc)
library(nlme)
library(corpcor)
library(doParallel)

source("combine_ute.R")
source("Exp_crit.R")
source("Laplace_approx.R")
source("log_posterior.R")
source("trace_mat.R")
source("Predict_uteApp.R")

cl <- makeCluster(9, outfile="cl_out.txt")
registerDoParallel(cl)

load("LP_post_New.RData")   # prior distribution
mu_prior <- LP_post[[1]]
Sigma_prior <- LP_post[[2]]   
iSigma_prior <- solve(Sigma_prior)
log_det_prior <- log(det(Sigma_prior))

N <- 5000 # number of samples from the prior
prior <- rmvnorm(N,mu_prior,Sigma_prior) # Initial particle set

num_par <- ncol(prior)  # total number of parameters in the Copula model
r01 <- 0.001            # nugget effect parameter for model 1
r02 <- 0.01             # nugget effect parameter for model 2

###################################################
Monitoring_sites <- read.csv("monitoring_sites_selected.csv") # prediction locations
Site_Locs <- data.frame(Monitoring_sites[,6:7]) 
dist1 <- rdist(Site_Locs)
dist_mat <- dist1/max(dist1)
tot_sites <- nrow(dist_mat)

data2017 <- read.csv("Meteorological_data2017.csv") # meteorological data at the prediction locations
X1 <- (-1*data2017$Loc_Y)/3000000
X2 <- data2017$Humidity/100
X3 <- data2017$Wind_speed/10
X.all <- cbind(X1,X2,X3)      # values of the covariates at the prediction locations


######## Optimal design selection using coordinate exchange algorithm ########

d_no <- 5   # number of design points
d_int <- sample(Monitoring_sites$Site_No,d_no,replace=F) # initial design

set.seed(101)
B1 <- 500
theta_int <- rmvnorm(n=B1,mean=mu_prior, sigma = Sigma_prior)

B2 <- 1000
Z1 <- matrix(rnorm(tot_sites*B2), tot_sites)

Unif_set <- runif(tot_sites,0,1)  
Z_res <- matrix(rnorm(tot_sites*3),ncol=3)

design_new <- list() 
ute_dist <- list()
ute_current <- (-10000)
d_prop <- d_int

for(k in 1:10){
  for(i in 1:d_no){
    Unsamp_new <- Monitoring_sites[!(Monitoring_sites$Site_No %in% d_prop[-i]),]
    ute <- foreach(j = 1:nrow(Unsamp_new),.packages= c("mvtnorm","corpcor","nlme","Matrix","matrixcalc","AirQualityRcpp","Rcpp","RcppArmadillo"),.verbose=TRUE,.combine = c) %dopar% 
    {
      d_new <- d_prop    
      d_new[i] <- Unsamp_new$Site_No[j]
      exp.crit(d=d_new,B=B1)
    }
    Ute_mat <- matrix(ute,ncol=nrow(Unsamp_new))
    max_d <- which.max(Ute_mat[1,])
    if(Ute_mat[1,max_d]>ute_current){
      d_prop[i] <- Unsamp_new$Site_No[max_d]
      ute_current <- Ute_mat[1,max_d]
      ute_dist[[(k-1)*d_no+i]] <- Ute_mat[,max_d]
    }else{
      ute_dist[[(k-1)*d_no+i]] <- ute_current
    }

  }
  
  design_new[[k]] <- d_prop

}  
  
out <- list(d_int,d_prop,ute_dist,design_new)

save(out, file = "out.RData")
stopCluster(cl)


