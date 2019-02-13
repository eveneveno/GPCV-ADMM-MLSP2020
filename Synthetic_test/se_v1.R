library(MASS)
library(tictoc)
library(crayon)
options(digits = 4)

# Kernel Function: SE kernel
myfun <- function(x,y){return((x-y)^2)}
K<- function(x,y,l){
  return(sig^2*exp(-outer(x,y,myfun)/(2*l^2)))
}

# Related Gradients
Kg <- function(x,y,l,sig=1){
  return(sig^2*exp(-outer(x,y,myfun)/(2*l^2)) * (outer(x,y,myfun)/l^3))
}

# Data generation
gen <- function(start,end,step,l,ratio,num_test){
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
  C = K(x,x,l)
  
  y = round(mvrnorm(1,rep(0,n),C) + runif(n,0,jit^2),4)
  yt = as.matrix(y[idx_train]); yv = as.matrix(y[idx_valid]); ys = as.matrix(y[idx_test])
  return(list(xt=xt,xv=xv,xs=xs,yt=yt,yv=yv,ys=ys))
}

# ADMM objective function
f <- function(l,z,lam){
  Cz = (K(xt,xt,l)+jit^2*diag(nt))%*%z
  return(norm(yv-K(xv,xt,l)%*%z,type="2")^2 + t(lam)%*%(Cz-yt) + rho/2*norm(Cz-yt,type="2")^2)
}

# Original objective function
partial1 <- function(l,z,lam){return(norm(yv-K(xv,xt,l)%*%z,type="2")^2)}

# Gradients
grad_l <- function(l,z,lam){
  Kg_tt = Kg(xt,xt,l); Kg_vt = Kg(xv,xt,l)
  K_tt = K(xt,xt,l); K_vt = K(xv,xt,l)
  
  grad = -2*t(yv)%*%Kg_vt%*%z + t(z)%*%t(Kg_vt)%*%K_vt%*%z + t(z)%*%t(K_vt)%*%Kg_vt%*%z + 
    t(lam)%*%Kg_tt%*%z + rho*t(jit^2*z-yt)%*%Kg_tt%*%z + (rho/2)*(t(z)%*%K_tt%*%Kg_tt%*%z + t(z)%*%Kg_tt%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_z <- function(l,z,lam){
  K_tt = K(xt,xt,l)
  K_vt = K(xv,xt,l)
  
  b <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt
  sigma <- 2*(crossprod(K_vt,K_vt%*%z) + 
                (rho/2)*(K_tt+jit^2*diag(length(xt))) %*% ((K_tt+jit^2*diag(length(xt)))%*%z))
  return(sigma + b)
}

#Plot result
ploty <- function(mc){
  layout(matrix(c(1,2,3,4,5,6,7,7,7,7,7,7),4,3,byrow = TRUE))
  
  #1.convergence of l
  plot(x=0:(length(l_collect)-1),y=l_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l')
  abline(h=l_true,col=3,lty=2)
  abline(h=l,col=4)
  text(round(l,5),x=(length(l_collect)-20), y=l+0.015)
  
  #2. auxilliary var
  plot(x=0:(length(z_collect)-1),y=z_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of z (l2 norm)')
  abline(h=norm(z,type="2"),col=4)
  text(round(norm(z,type="2"),6),x=(length(z_collect)-20), y=norm(z,type="2")+0.015) 
  
  #3. dual var
  plot(x=0:(length(lam_collect)-1),y=lam_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of lam (l2 norm)')
  abline(h=norm(lam,type="2"),col=4)
  text(round(norm(lam,type="2"),6),x=(length(lam_collect)-20), y=norm(lam,type="2")+0.015)
  
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
  
  # Fitting Performance
  plot(xt,yt,type='p',col='darkgrey',cex=0.5,xlab='x',ylab='fitted value')
  points(xv,yv,cex=0.5)
  
  y_pred = K(xv,xt,l)%*%z #(1) Convergent
  lines(xv,y_pred,col=2)
  
  z_true = solve(K(xt,xt,l_true)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,l_true)%*%z_true #(2) TRUE
  lines(xv,y_pred,col=3)
  
  z_init = solve(K(xt,xt,l_init)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,l_init)%*%z_init #(3) Initial
  lines(xv,y_pred,col=5,lty=2)
  
  legend(0,max(c(yt,yv)), legend=c("Prediction with the underlying truth", "Prediction with convergent params","Predictions with initial params"),
         col=c(3,2,5),lty=1, cex=0.8)
  legend(0,min(c(yt,yv))+1, legend=c("training pts", "validation pts"),
         col=c('darkgrey','black'),pch=c(1,1), cex=0.8)
}

#kernel setting
l_true  = 0.5

jit = 0.5 
sig = 1

alpha = beta = 0.5
rho = 5
rhorho = 5
Maxiter = 10

Num_mc = 10
l_corpora = c()

for (mc in 1:Num_mc){
  data = gen(0,10,0.02,l_true,1/2,20)
  xt = data$xt; xv = data$xv; xs = data$xs; yt = data$yt; yv = data$yv; ys = data$ys; nt = length(xt)
  # file = paste(format(Sys.time(), "%b_%d %X"),mc)
  # write.csv(data.frame(xt=xt,yt=yt),file=paste(file,'train.csv'),row.names = F)
  # write.csv(data.frame(xv=xv,yv=yv),file=paste(file,'valid.csv'),row.names = F)
  # write.csv(data.frame(xs=xs,ys=ys),file=paste(file,'test.csv'),row.names = F)
  
  ##====================1-fold=======================##
  l_collect = c();cost_collect = c();part1_collect = c()

  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l = l_init = 0.9
  l_collect = l
  
  l1 = abs(round(rnorm(1,mean = l_init,sd = 0.1),3))
  z = z_init = solve(K(xt,xt,l1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l)+jit^2*diag(nt),yt),type="2")
  f_collect = c()
  p1_collect = c()
  
  cat("MC:",mc,"RS:",k," | l init:",green(l),"| jitter",green(l1),
      blue(round(f(l,z,lam),4)),blue(round(partial1(l,z,lam),4)),'\n')
  
  for (i in 1:Maxiter){
    tic()
    ##==l1 update==##
    iter = 0
    while(norm(grad_l(l,z,lam),type="2")>1e-03){
      
      d = grad_l(l,z,lam)
      t = 0.125
      while(f(l-t*d,z,lam)>=f(l,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l = l-t*d
      iter = iter+1
    }
    l_collect = c(l_collect,l)
    f_collect = c(f_collect,f(l,z,lam))
    p1_collect = c(p1_collect,partial1(l,z,lam))
    cat(i,red(round(l,4)),yellow(iter),blue(round(f(l,z,lam),3)),'|')
    
    ##==Z update==##
    iter = length(z)
    K_tt = K(xt,xt,l)
    K_vt = K(xv,xt,l)
    bv <- -2*crossprod(K_vt,yv) + 
      (K_tt+jit^2*diag(length(xt))) %*% lam - 
      rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
    sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
    
    d = -grad_z(l,z,lam)
    for (ii in 1:(iter-1)){
      stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
      zold = z
      z = z+stepsize*d
      if(norm(z-zold,type="2") <= 1e-5){break}
      bb = as.numeric((crossprod(grad_z(l,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
      d = -grad_z(l,z,lam)+bb*d
    }
    
    z_collect = c(z_collect,norm(z,type='2'))
    f_collect = c(f_collect,f(l,z,lam))
    p1_collect = c(p1_collect,partial1(l,z,lam))
    cat(norm(z,type='2'),yellow(ii),blue(round(f(l,z,lam),5)))
    
    gap = norm(z-solve(K(xt,xt,l)+jit^2*diag(nt),yt),type="2")
    gap_collect = c(gap_collect,gap)
    lam = lam + rhorho*((K(xt,xt,l)+jit^2*diag(nt))%*%z - yt)
    cat('GAP',gap,' |',blue(round(f(l,z,lam),4)),' |')
    lam_collect = c(lam_collect,norm(lam,type='2'))
    f_collect = c(f_collect,f(l,z,lam))
    p1_collect = c(p1_collect,partial1(l,z,lam))
    
    toc()
  }
  l_left = l
  ploty(mc)
  
  
  ##====================2-fold=======================##
  
  xt = data$xv; xv = data$xt; yv = data$yt; yt = data$yv; nt = length(xt)
  l_collect = c();cost_collect = c();part1_collect = c()
  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l = l_init = 0.9
  l_collect = l
  
  l1 = abs(round(rnorm(1,mean = l_init,sd = 0.1),3))
  z = z_init = solve(K(xt,xt,l1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l)+jit^2*diag(nt),yt),type="2")
  f_collect = c()
  p1_collect = c()
  
  cat("MC:",mc,"RS:",k," | l init:",green(l),"| jitter",green(l1),
      blue(round(f(l,z,lam),4)),blue(round(partial1(l,z,lam),4)),'\n')
  
  for (k in 1:Num_restart){
    lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
    l = l_init = 1.5
    l1 = abs(round(rnorm(1,mean = l_init,sd = 0.1),3))
    l_collect = l
    z = z_init = solve(K(xt,xt,l1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
    
    gap_collect = norm(z-solve(K(xt,xt,l)+jit^2*diag(nt),yt),type="2")
    f_collect = c()
    p1_collect = c()
    
    cat("MC:",mc,"RS:",k," | l init:",green(l),"| jitter",green(l),
        blue(round(f(l,z,lam),4)),blue(round(partial1(l,z,lam),4)),'\n')
    
    for (i in 1:Maxiter){
      tic()
      ##==l1 update==##
      iter = 0
      while(norm(grad_l(l,z,lam),type="2")>1e-03){
        
        d = grad_l(l,z,lam)
        t = 0.125
        while(f(l-t*d,z,lam)>=f(l,z,lam)-alpha*t*d*d){
          t = t*beta
        }
        l = l-t*d
        iter = iter+1
      }
      l_collect = c(l_collect,l)
      f_collect = c(f_collect,f(l,z,lam))
      p1_collect = c(p1_collect,partial1(l,z,lam))
      cat(i,red(round(l,4)),yellow(iter),blue(round(f(l,z,lam),3)),'|')
      
      
      ##==Z update==##
      iter = length(z)
      K_tt = K(xt,xt,l)
      K_vt = K(xv,xt,l)
      bv <- -2*crossprod(K_vt,yv) + 
        (K_tt+jit^2*diag(length(xt))) %*% lam - 
        rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
      sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
      
      d = -grad_z(l,z,lam)
      for (ii in 1:(iter-1)){
        stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
        zold = z
        z = z+stepsize*d
        if(norm(z-zold,type="2") <= 1e-5){break}
        bb = as.numeric((crossprod(grad_z(l,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
        d = -grad_z(l,z,lam)+bb*d
      }
      
      z_collect = c(z_collect,norm(z,type='2'))
      f_collect = c(f_collect,f(l,z,lam))
      p1_collect = c(p1_collect,partial1(l,z,lam))
      cat(norm(z,type='2'),yellow(ii),blue(round(f(l,z,lam),5)))
      
      gap = norm(z-solve(K(xt,xt,l)+jit^2*diag(nt),yt),type="2")
      gap_collect = c(gap_collect,gap)
      lam = lam + rhorho*((K(xt,xt,l)+jit^2*diag(nt))%*%z - yt)
      cat('GAP',gap,' |',blue(round(f(l,z,lam),4)),' |')
      lam_collect = c(lam_collect,norm(lam,type='2'))
      f_collect = c(f_collect,f(l,z,lam))
      p1_collect = c(p1_collect,partial1(l,z,lam))
      
      toc()
    }
    l_right = l
    ploty(mc,k)
  }
  ##===============average===================##
  l_right = l1
  l_corpora = c(l_corpora,(l_left+l_right)/2)
}
# file = paste(format(Sys.time(), "%b_%d %X"),mc)
# write.csv(data.frame(l=l_corpora,mse1=mse1_corpora,mse2=mse2_corpora,mse=mse_corpora,mse_truth=mse_t_corpora),file=paste(file,'_params.csv',sep=''),row.names = F)
