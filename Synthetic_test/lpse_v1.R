library(MASS)
library(tictoc)
library(crayon)
options(digits = 4)

# Kernel Function: LP+SE kernel
myfun_per <- function(x,y) {return(pi*abs(x-y))}
myfun_se <- function(x,y) {return((x-y)^2)}

K<- function(x,y,l1,l2,p){
  inn_per = outer(x,y,myfun_per)
  inn_se = outer(x,y,myfun_se)
  return(exp(-inn_se/(2*l1^2)) + exp(-2*(sin(inn_per/p)^2)/(l2^2))*exp(-inn_se/(2*l2^2)))
}

Kg_l1 <- function(x,y,l1,l2,p){
  inn_se = outer(x,y,myfun_se)
  return(exp(-inn_se/(2*l1^2)) * (inn_se/l1^3))
}

Kg_l2 <- function(x,y,l1,l2,p){
  inn_se = outer(x,y,myfun_se)
  inn_per = outer(x,y,myfun_per)
  return(exp(-inn_se/(2*l2^2))*(inn_se/l2^3) * exp(-2*(sin(inn_per/p)^2)/(l2^2)) + 
           exp(-2*(sin(inn_per/p)^2)/(l2^2))*(4*(sin(inn_per/p)^2)/(l2^3))*exp(-inn_se/(2*l2^2)))
}

Kg_p <- function(x,y,l1,l2,p){
  inn_se = outer(x,y,myfun_se)
  inn_per = outer(x,y,myfun_per)
  return(exp(-2*(sin(inn_per/p)^2)/(l2^2)) * (-4*sin(inn_per/p)/(l2^2)) * cos(inn_per/p) * (-inn_per/(p^2))*exp(-inn_se/(2*l2^2)))
}

# Data generation
gen <- function(start,end,step,l1,l2,p,ratio,num_test){
  x = seq(start,end,step)
  n = length(x)
  
  idx_train = seq(1,n,2) 
  idx_test_t = sort(sample(idx_train,num_test/2))
  idx_train = idx_train[-idx_test_t]
  
  idx_valid = seq(2,n,2)
  idx_test_v = sort(sample(idx_valid,num_test/2))
  idx_valid = idx_valid[-idx_test_v]
  
  idx_test = sort(c(idx_test_t,idx_test_v))
  
  xt = x[idx_train]; xv = x[idx_valid]; xs = x[idx_test]
  C = K(x,x,l1,l2,p)
  
  y = round(mvrnorm(1,rep(0,n),C) + runif(n,0,jit^2),4)
  yt = as.matrix(y[idx_train]); yv = as.matrix(y[idx_valid]); ys = as.matrix(y[idx_test])
  return(list(xt=xt,xv=xv,xs=xs,yt=yt,yv=yv,ys=ys))
}

# ADMM objective function
f <- function(l1,l2,p,z,lam){
  Cz = (K(xt,xt,l1,l2,p)+jit^2*diag(nt))%*%z
  return(norm(yv-K(xv,xt,l1,l2,p)%*%z,type="2")^2 + t(lam)%*%(Cz-yt) + rho/2*norm(Cz-yt,type="2")^2)
}

# Original objective function
partial1 <- function(l1,l2,p,z,lam){return(norm(yv-K(xv,xt,l1,l2,p)%*%z,type="2")^2)}

# Gradients
grad_l1 <- function(l1,l2,p,z,lam){
  Kg_tt = Kg_l1(xt,xt,l1,l2,p); Kg_vt = Kg_l1(xv,xt,l1,l2,p)
  K_tt = K(xt,xt,l1,l2,p); K_vt = K(xv,xt,l1,l2,p)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z +crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_l2 <- function(l1,l2,p,z,lam){
  Kg_tt = Kg_l2(xt,xt,l1,l2,p); Kg_vt = Kg_l2(xv,xt,l1,l2,p)
  K_tt = K(xt,xt,l1,l2,p); K_vt = K(xv,xt,l1,l2,p)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z + crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_p <- function(l1,l2,p,z,lam){
  Kg_tt = Kg_p(xt,xt,l1,l2,p); Kg_vt = Kg_p(xv,xt,l1,l2,p)
  K_tt = K(xt,xt,l1,l2,p); K_vt = K(xv,xt,l1,l2,p)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z + crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_z <- function(l1,l2,p,z,lam){
  K_tt = K(xt,xt,l1,l2,p)
  K_vt = K(xv,xt,l1,l2,p)
  
  b <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt
  sigma <- 2*(crossprod(K_vt,K_vt%*%z) + 
                (rho/2)*(K_tt+jit^2*diag(length(xt))) %*% ((K_tt+jit^2*diag(length(xt)))%*%z))
  return(sigma + b)
}

#Plot result
ploty <- function(mc,k){
  layout(matrix(c(1,2,3,4,5,6,7,7,7,7,7,7),4,3,byrow = TRUE))
  #1.
  plot(x=0:(length(l1_collect)-1),y=l1_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l1')
  abline(h=l1_true,col=3,lty=2)
  abline(h=l1,col=4)
  text(round(l1,5),x=(length(l1_collect)-20), y=l1+0.015)
  
  #2.
  plot(x=0:(length(l2_collect)-1),y=l2_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l2')
  abline(h=l2_true,col=3,lty=2)
  abline(h=l2,col=4)
  text(round(l2,5),x=(length(l2_collect)-20), y=l2+0.015)
  
  #3.
  plot(x=0:(length(p_collect)-1),y=p_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of p')
  abline(h=p_true,col=3,lty=2)
  abline(h=p,col=4)
  text(round(p,5),x=(length(p_collect)-20), y=p+0.015)
  
  #4.Equality Gap
  plot(x=0:(length(gap_collect)-1),y=gap_collect,type='l',col=2,
       xlab = 'Iter',ylab='Inequality Gap Control')
  text(paste('Rho:',rho),x=(length(gap_collect)/2), y=(max(gap_collect)+min(gap_collect))/1.8) 
  text(paste('RhoRho:',rhorho),x=(length(gap_collect)/2), y=(max(gap_collect)+min(gap_collect))/2.2) 
  
  #5.ADMM obj
  plot(f_collect,type='l',cex=0.5,col=2,ylab='ADMM objective')
  abline(h=tail(f_collect,n=1),col=4)
  text(round(tail(f_collect,n=1),6),x=(length(f_collect)-20), y=tail(f_collect, n=1)+0.015) 
  
  #6.Original obj
  plot(p1_collect,type='l',cex=0.5,col=2,ylab='Original objective')
  abline(h=tail(p1_collect, n=1),col=4)    ###==== When rho is small, may display unusual problem ====###
  text(round(tail(p1_collect,n=1),6),x=(length(p1_collect)-20), y=tail(p1_collect, n=1)+0.015) 
  
  #4.fitted value
  plot(xt,yt,type='p',col='darkgrey',cex=0.5,xlab='x',ylab='fitted value')
  points(xv,yv,cex=0.5)
  
  y_pred = K(xv,xt,l1,l2,p)%*%z #(1) Convergent
  lines(xv,y_pred,col=2)
  
  z_true = solve(K(xt,xt,l1_true,l2_true,p_true)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,l1_true,l2_true,p_true)%*%z_true #(2) TRUE
  lines(xv,y_pred,col=3)
  
  z_init = solve(K(xt,xt,l1_init,l2_init,p_init)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,l1_init,l2_true,p_init)%*%z_init #(3) Initial
  lines(xv,y_pred,col=5,lty=2)
  
  legend(0,max(c(yt,yv)), legend=c("Prediction with the underlying truth", "Prediction with convergent params","Predictions with initial params"),
         col=c(3,2,5),lty=1, cex=0.8)
  legend(0,min(c(yt,yv))+1, legend=c("training pts", "validation pts"),
         col=c('darkgrey','black'),pch=c(1,1), cex=0.8)
}

# kernel setting
l1_true = 3
l2_true = 1
p_true = 2

jit = 0.5 
sig = 1

alpha = beta = 0.5
rho = 5
rhorho = 5
Maxiter = 100

Num_mc = 10
l1_corpora = c()
l2_corpora = c()
p_corpora = c()

for (mc in 1:Num_mc){
  data = gen(0,10,0.03,l1_true,l2_true,p_true,1/2,20)
  xt = data$xt; xv = data$xv; xs = data$xs; yt = data$yt; yv = data$yv; ys = data$ys; nt = length(xt)
  #file = paste(format(Sys.time(), "%b_%d %X"),mc)
  #write.csv(data.frame(xt=xt,yt=yt),file=paste(file,'train.csv'),row.names = F)
  #write.csv(data.frame(xv=xv,yv=yv),file=paste(file,'valid.csv'),row.names = F)
  #write.csv(data.frame(xs=xs,ys=ys),file=paste(file,'test.csv'),row.names = F)
  
  #1-fold
  l1_collect = c(); l2_collect = c(); p_collect = c()
  cost_collect = c();part1_collect = c() 
  
  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l1 = l1_init = 3.5 
  l2 = l2_init = 1.5 
  p = p_init = 2.5 
  l11 = abs(round(rnorm(1,mean = l1_init,sd = 0.1),3))
  l21 = abs(round(rnorm(1,mean = l2_init,sd = 0.1),3))
  p1 = abs(round(rnorm(1,mean = p_init,sd = 0.1),3))
  l1_collect = l1; l2_collect = l2; p_collect = p
  z = z_init = solve(K(xt,xt,l1,l2,p1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l1,l2,p)+jit^2*diag(nt),yt),type="2")
  f_collect = c()
  p1_collect = c()
  
  cat("MC:",mc,"RS:",k," | l1 l2 p init:",green(l1,l2,p),"| jitter",green(l11,l21,p1),
      blue(round(f(l1,l2,p,z,lam),4)),blue(round(partial1(l1,l2,p,z,lam),4)),'\n')
    
  for (i in 1:Maxiter){
    tic()
    ##==l1 update==##
    iter = 0
    while(norm(grad_l1(l1,l2,p,z,lam),type="2")>1e-03){
      
      d = grad_l1(l1,l2,p,z,lam)
      t = 0.125
      while(f(l1-t*d,l2,p,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l1 = l1-t*d
      iter = iter+1
    }
    l1_collect = c(l1_collect,l1)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(i,red(round(l1,4)),yellow(iter),blue(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==l2 update==##
    iter = 0
    while(norm(grad_l2(l1,l2,p,z,lam),type="2")>1e-03){
      d = grad_l2(l1,l2,p,z,lam)
      t = 0.0625
      while(f(l1,l2-t*d,p,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l2 = l2-t*d
      iter = iter+1
    }
    l2_collect = c(l2_collect,l2)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(red(round(l2,4)),yellow(iter),blue(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==p update==##
    iter = 0
    while(norm(grad_p(l1,l2,p,z,lam),type="2")>1e-03){
      d = grad_p(l1,l2,p,z,lam)
      t = 0.0625
      while(f(l1,l2,p-t*d,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      p = p-t*d
      iter = iter+1
    }
    p_collect = c(p_collect,p)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(red(round(p,4)),yellow(iter),cyan(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==Z update==##
    iter = length(z)
    K_tt = K(xt,xt,l1,l2,p)
    K_vt = K(xv,xt,l1,l2,p)
    bv <- -2*crossprod(K_vt,yv) + 
      (K_tt+jit^2*diag(length(xt))) %*% lam - 
      rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
    sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
    
    d = -grad_z(l1,l2,p,z,lam)
    for (ii in 1:(iter-1)){
      stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
      zold = z
      z = z+stepsize*d
      if(norm(z-zold,type="2") <= 1e-5){break}
      bb = as.numeric((crossprod(grad_z(l1,l2,p,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
      d = -grad_z(l1,l2,p,z,lam)+bb*d
    }
    
    z_collect = c(z_collect,norm(z,type='2'))
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(norm(z,type='2'),yellow(ii),blue(round(f(l1,l2,p,z,lam),5)))
    
    gap = norm(z-solve(K(xt,xt,l1,l2,p)+jit^2*diag(nt),yt),type="2")
    gap_collect = c(gap_collect,gap)
    lam = lam + rhorho*((K(xt,xt,l1,l2,p)+jit^2*diag(nt))%*%z - yt)
    cat('GAP',gap,' |',blue(round(f(l1,l2,p,z,lam),4)),' |')
    lam_collect = c(lam_collect,norm(lam,type='2'))
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    toc()
  }
  l1_left = l1; l2_left = l2; p_left = p
  ploty(mc)
  
  #2-fold
  xt = data$xv; xv = data$xt; yt = data$yv; yv = data$yt;nt = length(xt)
  l1_collect = c(); l2_collect = c(); p_collect = c()
  cost_collect = c();part1_collect = c() 
  
  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l1 = l1_init = 3.5 
  l2 = l2_init = 1.5 
  p = p_init = 2.5 
  l11 = abs(round(rnorm(1,mean = l1_init,sd = 0.1),3))
  l21 = abs(round(rnorm(1,mean = l2_init,sd = 0.1),3))
  p1 = abs(round(rnorm(1,mean = p_init,sd = 0.1),3))
  l1_collect = l1; l2_collect = l2; p_collect = p
  z = z_init = solve(K(xt,xt,l1,l2,p1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l1,l2,p)+jit^2*diag(nt),yt),type="2")
  f_collect = c()
  p1_collect = c()
  
  cat("MC:",mc,"RS:",k," | l1 l2 p init:",green(l1,l2,p),"| jitter",green(l11,l21,p1),
      blue(round(f(l1,l2,p,z,lam),4)),blue(round(partial1(l1,l2,p,z,lam),4)),'\n')
  
  for (i in 1:Maxiter){
    tic()
    ##==l1 update==##
    iter = 0
    while(norm(grad_l1(l1,l2,p,z,lam),type="2")>1e-03){
      
      d = grad_l1(l1,l2,p,z,lam)
      t = 0.125
      while(f(l1-t*d,l2,p,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l1 = l1-t*d
      iter = iter+1
    }
    l1_collect = c(l1_collect,l1)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(i,red(round(l1,4)),yellow(iter),blue(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==l2 update==##
    iter = 0
    while(norm(grad_l2(l1,l2,p,z,lam),type="2")>1e-03){
      d = grad_l2(l1,l2,p,z,lam)
      t = 0.0625
      while(f(l1,l2-t*d,p,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l2 = l2-t*d
      iter = iter+1
    }
    l2_collect = c(l2_collect,l2)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(red(round(l2,4)),yellow(iter),blue(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==p update==##
    iter = 0
    while(norm(grad_p(l1,l2,p,z,lam),type="2")>1e-03){
      d = grad_p(l1,l2,p,z,lam)
      t = 0.0625
      while(f(l1,l2,p-t*d,z,lam)>=f(l1,l2,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      p = p-t*d
      iter = iter+1
    }
    p_collect = c(p_collect,p)
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(red(round(p,4)),yellow(iter),cyan(round(f(l1,l2,p,z,lam),3)),'|')
    
    ##==Z update==##
    iter = length(z)
    K_tt = K(xt,xt,l1,l2,p)
    K_vt = K(xv,xt,l1,l2,p)
    bv <- -2*crossprod(K_vt,yv) + 
      (K_tt+jit^2*diag(length(xt))) %*% lam - 
      rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
    sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
    
    d = -grad_z(l1,l2,p,z,lam)
    for (ii in 1:(iter-1)){
      stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
      zold = z
      z = z+stepsize*d
      if(norm(z-zold,type="2") <= 1e-5){break}
      bb = as.numeric((crossprod(grad_z(l1,l2,p,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
      d = -grad_z(l1,l2,p,z,lam)+bb*d
    }
    
    z_collect = c(z_collect,norm(z,type='2'))
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    cat(norm(z,type='2'),yellow(ii),blue(round(f(l1,l2,p,z,lam),5)))
    
    gap = norm(z-solve(K(xt,xt,l1,l2,p)+jit^2*diag(nt),yt),type="2")
    gap_collect = c(gap_collect,gap)
    lam = lam + rhorho*((K(xt,xt,l1,l2,p)+jit^2*diag(nt))%*%z - yt)
    cat('GAP',gap,' |',blue(round(f(l1,l2,p,z,lam),4)),' |')
    lam_collect = c(lam_collect,norm(lam,type='2'))
    f_collect = c(f_collect,f(l1,l2,p,z,lam))
    p1_collect = c(p1_collect,partial1(l1,l2,p,z,lam))
    toc()
  }
  l1_right = l1; l2_right = l2; p_right = p
  ploty(mc)
  
  l1_corpora = c(l1_corpora,(l1_left+l1_right)/2)
  l2_corpora = c(l2_corpora,(l2_left+l2_right)/2)
  p_corpora = c(p_corpora,(p_left+p_right)/2)
}

#file = paste(format(Sys.time(), "%b_%d %X"),mc)
#write.csv(data.frame(l=l_corpora,p=p_corpora),file=paste(file,'_params.csv',sep=''),row.names = F)