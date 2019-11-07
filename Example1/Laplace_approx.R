Laplace_approx<- function(X,Y,theta,distM)
  {
  #lp.approx <- optim(par=theta,X=X, Y=Y,distM=distM, fn=log_post, hessian=TRUE)
  lp.approx <-nlminb(start=theta,X=X,distM=distM, Y=Y, objective =log_post,lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max))) 
  
  lp.approx
}