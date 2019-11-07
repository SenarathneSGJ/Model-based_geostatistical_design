exp.crit <- function(d,B)
{
  distM=rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)

  crit=foreach(j = 1:B,.packages= c("mvtnorm","fields","corpcor","nlme","Matrix","matrixcalc","SpatialSimEx","RcppArmadillo"),.errorhandling = 'stop',.export= ls(globalenv()),.combine = c) %dopar% 
  {
    Y <- Res_Clcpp(xdata=as.matrix(X),distM,para=theta_int[j,],Z1=t(Z_res[1:nrow(X),1]),Z2=t(Z_res[1:nrow(X),2]),uv=Unif_set[1:nrow(X),],r01,r02)
    crit.out <- combine_ute(d,X,Y,theta=theta_int[j,],distM)
    crit_both <- crit.out[[3]]      # dual purpose utility
    #crit_kld <- crit.out[[4]]      # estimation utility
    #crit_pred <- crit.out[[5]]     # prediction utility
    crit_both
  }
  
  return(crit)
}
