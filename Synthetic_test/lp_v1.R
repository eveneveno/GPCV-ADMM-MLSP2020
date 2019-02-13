library(MASS)
library(tictoc)
library(crayon)
options(digits = 4)

# Kernel Function: LP kernel
myfun_per <- function(x,y) {return(pi*abs(x-y))}
myfun_se <- function(x,y) {return((x-y)^2)}

K<- function(x,y,l,p){
  inn_per = outer(x,y,myfun_per);inn_se = outer(x,y,myfun_se)
  return(sig^2*exp(-2*(sin(inn_per/p)^2)/(l^2))*exp(-inn_se/(2*l^2)))
}

# Related Gradients
Kg_l <- function(x,y,l,p){
  inn_se = outer(x,y,myfun_se)
  inn_per = outer(x,y,myfun_per)
  return(sig^2*exp(-inn_se/(2*l^2))*(inn_se/l^3) * exp(-2*(sin(inn_per/p)^2)/(l^2)) + 
           sig^2*exp(-2*(sin(inn_per/p)^2)/(l^2))*(4*(sin(inn_per/p)^2)/(l^3))*exp(-inn_se/(2*l^2)))
}

Kg_p <- function(x,y,l,p){
  inn_se = outer(x,y,myfun_se)
  inn_per = outer(x,y,myfun_per)
  return(sig^2*exp(-2*(sin(inn_per/p)^2)/(l^2)) * (-4*sin(inn_per/p)/(l^2)) * cos(inn_per/p) * (-inn_per/(p^2))*exp(-inn_se/(2*l^2)))
}

# Data generation
gen <- function(start,end,step,l,p,ratio,num_test){
  x = seq(start,end,step)
  n = length(x)
  
  #let odd numbers be training points
  idx_train = seq(1,n,2) 
  idx_test_t = sort(sample(idx_train,num_test/2)) #separate half of test data from odd numbers
  idx_train = idx_train[-idx_test_t]
  
  #let even numbers be validation points
  idx_valid = seq(2,n,2) 
  idx_test_v = sort(sample(idx_valid,num_test/2)) #separate half of test data from even numbers
  idx_valid = idx_valid[-idx_test_v]
  
  idx_test = sort(c(idx_test_t,idx_test_v))
  
  xt = x[idx_train]; xv = x[idx_valid]; xs = x[idx_test]
  C = K(x,x,l,p)
  
  y = round(mvrnorm(1,rep(0,n),C) + runif(n,0,jit^2),4) #generate y values
  yt = as.matrix(y[idx_train]); yv = as.matrix(y[idx_valid]); ys = as.matrix(y[idx_test])
  return(list(xt=xt,xv=xv,xs=xs,yt=yt,yv=yv,ys=ys))
}

# ADMM objective function
f <- function(l,p,z,lam){
  Cz = (K(xt,xt,l,p)+jit^2*diag(nt))%*%z
  return(norm(yv-K(xv,xt,l,p)%*%z,type="2")^2 + t(lam)%*%(Cz-yt) + rho/2*norm(Cz-yt,type="2")^2)
}

# Original objective function
partial1 <- function(l,p,z,lam){return(norm(yv-K(xv,xt,l,p)%*%z,type="2")^2)}

# Gradients
grad_l <- function(l,p,z,lam){
  Kg_tt = Kg_l(xt,xt,l,p); Kg_vt = Kg_l(xv,xt,l,p)
  K_tt = K(xt,xt,l,p); K_vt = K(xv,xt,l,p)
  
  grad = -2*t(yv)%*%Kg_vt%*%z + t(z)%*%t(Kg_vt)%*%K_vt%*%z + t(z)%*%t(K_vt)%*%Kg_vt%*%z + 
    t(lam)%*%Kg_tt%*%z + rho*t(jit^2*z-yt)%*%Kg_tt%*%z + (rho/2)*(t(z)%*%K_tt%*%Kg_tt%*%z + t(z)%*%Kg_tt%*%K_tt%*%z)
  return(as.numeric(grad))
}

#gradient of ADMM objective function wrt p
grad_p <- function(l,p,z,lam){
  Kg_tt = Kg_p(xt,xt,l,p); Kg_vt = Kg_p(xv,xt,l,p)
  K_tt = K(xt,xt,l,p); K_vt = K(xv,xt,l,p)
  
  grad = -2*t(yv)%*%Kg_vt%*%z + t(z)%*%t(Kg_vt)%*%K_vt%*%z + t(z)%*%t(K_vt)%*%Kg_vt%*%z + 
    t(lam)%*%Kg_tt%*%z + rho*t(jit^2*z-yt)%*%Kg_tt%*%z + (rho/2)*(t(z)%*%K_tt%*%Kg_tt%*%z + t(z)%*%Kg_tt%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_z <- function(l,p,z,lam){
  K_tt = K(xt,xt,l,p)
  K_vt = K(xv,xt,l,p)
  
  b <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt
  sigma <- 2*(crossprod(K_vt,K_vt%*%z) + 
                (rho/2)*(K_tt+jit^2*diag(length(xt))) %*% ((K_tt+jit^2*diag(length(xt)))%*%z))
  return(sigma + b)
}

#Plot result
ploty <- function(mc){
  
  par(mar=c(2,2,2,1))
  layout(matrix(c(1,2,3,4,5),1,5,byrow = TRUE))
  
  #1.convergence of l
  plot(x=0:(length(l_collect)-1),y=l_collect,type='l',col=2,
       ylab='',xlab='',lwd=6,cex.lab=1.5,font.lab=2)
  abline(h=l,col=4,lwd=4)
  title('Hyper-parameter l',font=2)
  
  #2.convergence of p
  plot(x=0:(length(p_collect)-1),y=p_collect,type='l',col=2,
       ylab='',xlab='',lwd=6,cex.lab=1.5,font.lab=2)
  abline(h=p,col=4,lwd=4)
  title('Hyper-parameter p',font=2)  

  #3.auxilliary variable norm
  plot(x=0:(length(z_collect)-1),y=z_collect,type='l',col='salmon',lwd=6,
       ylab='',xlab='',cex.lab=1.5,font.lab=2)
  abline(h=norm(z,type="2"),col=4,lwd=4)
  title('Auxiliary Variable Norm',font=2)
  
  # dual var
  # plot(x=0:(length(lam_collect)-1),y=lam_collect,type='l',col='turquoise',lwd=6,
  #      xlab = 'Iter',ylab='change of lam (l2 norm)')
  # abline(h=norm(lam,type="2"),col=4,lwd=4)
  # text(round(norm(lam,type="2"),6),x=(length(lam_collect)-20), y=norm(lam,type="2")+0.015)
  # title('Dual Variable Norm')
  
  #4.Equality Gap
  plot(x=0:(length(gap_collect)-1),y=gap_collect,type='l',col='turquoise',
       ylab='',xlab='',lwd=6,cex.lab=1.5,font.lab=2)
  abline(h=tail(gap_collect,n=1),col=4,lwd=4)
  title('Inequality Gap',font=2)
  
  #5.Original obj
  plot(p1_collect,type='l',cex=0.5,col='slateblue',ylab='Objective value',xlab='',lwd=6,cex.lab=1.5,font.lab=2)
  abline(h=tail(p1_collect, n=1),col=4,lwd=4)  
  title('Objective Value',font=2)
  
  par(mfrow=c(1,1))
  # Fitting Performance
  plot(xt,yt,type='p',cex=1.8,pch=16,xlab='',ylab='')
  points(xv,yv,cex=2,pch=16)
  
  y_pred = K(xv,xt,l,p)%*%z #(1) Convergent
  lines(xv,y_pred,col=2,lwd=8)
  
  z_true = solve(K(xt,xt,l_true,p_true)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,l_true,p_true)%*%z_true #(2) TRUE
  lines(xv,y_pred,col=3,lwd=4)
  
  z_init = solve(K(xt,xt,1,2)+jit^2*diag(nt),yt)
  y_pred = K(xv,xt,1,2)%*%z_init #(3) Initial
  lines(xv,y_pred,col=5,lty=4,lwd=5)
} 

#kernel setting
l_true = 0.5
p_true = 1

jit = 0.5 #noise term
sig = 1 #sigma

alpha = beta = 0.5
rho = 5
rhorho = 5
Maxiter = 20 

Num_mc = 1
l_corpora = c()
p_corpora = c()

for (mc in 1:Num_mc){
  data = gen(0,20.4,0.02,l_true,p_true,1/2,20) #500 for training, 20 for testing
  xt = data$xt; xv = data$xv; xs = data$xs; yt = data$yt; yv = data$yv; ys = data$ys; nt = length(xt)
  #file = paste(format(Sys.time(), "%b_%d %X"),mc)
  #write.csv(data.frame(xt=xt,yt=yt),file=paste(file,'train.csv'),row.names = F)
  #write.csv(data.frame(xv=xv,yv=yv),file=paste(file,'valid.csv'),row.names = F)
  #write.csv(data.frame(xs=xs,ys=ys),file=paste(file,'test.csv'),row.names = F)
  
  #1-fold
  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l = l_init = 0.9 
  p = p_init = 1.5
  l_collect = l
  p_collect = p
  
  l1 = abs(round(rnorm(1,mean = l_init,sd = 0.1),3))
  p1 = abs(round(rnorm(1,mean = p_init,sd = 0.1),3))
  z = z_init = solve(K(xt,xt,l1,p1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l,p)+jit^2*diag(nt),yt),type="2")
  
  f_collect = c();p1_collect = c()
  
  cat("MC:",mc,"| l init:",green(l,p),"| jitter",green(l1,p1),
      blue(round(f(l,p,z,lam),4)),blue(round(partial1(l,p,z,lam),4)),'\n')
  
  for (i in 1:Maxiter){
    tic()
    ##==l update==##
    iter = 0
    while(norm(grad_l(l,p,z,lam),type="2")>1e-03){
      d = grad_l(l,p,z,lam)
      t = 0.125
      while(f(l-t*d,p,z,lam)>=f(l,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l = l-t*d
      iter = iter+1
    }
    l_collect = c(l_collect,l)
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat(i,'l:',red(round(l,4)),yellow(iter),blue(round(f(l,p,z,lam),3)),'|')
    
    ##==p update==##
    iter = 0
    while(norm(grad_p(l,p,z,lam),type="2")>1e-03){
      d = grad_p(l,p,z,lam)
      t = 0.125
      while(f(l,p-t*d,z,lam)>=f(l,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      p = p-t*d
      iter = iter+1
    }
    p_collect = c(p_collect,p)
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat('p:',red(round(p,4)),yellow(iter),blue(round(f(l,p,z,lam),3)),'|')
    
    ##==Z update==##
    iter = length(z)
    K_tt = K(xt,xt,l,p)
    K_vt = K(xv,xt,l,p)
    bv <- -2*crossprod(K_vt,yv) + 
      (K_tt+jit^2*diag(length(xt))) %*% lam - 
      rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
    sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
    
    d = -grad_z(l,p,z,lam)
    for (ii in 1:(iter-1)){
      stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
      zold = z
      z = z+stepsize*d
      if(norm(z-zold,type="2") <= 1e-5){break}
      bb = as.numeric((crossprod(grad_z(l,p,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
      d = -grad_z(l,p,z,lam)+bb*d
    }
    
    z_collect = c(z_collect,norm(z,type='2'))
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat('z norm:',red(norm(z,type='2')),yellow(ii),blue(round(f(l,p,z,lam),5)))
    
    gap = norm(z-solve(K(xt,xt,l,p)+jit^2*diag(nt),yt),type="2");gap_collect = c(gap_collect,gap)
    cat('| GAP',gap,'|',blue(round(f(l,p,z,lam),4)),'|')
    
    #dual update
    lam = lam + rhorho*((K(xt,xt,l,p)+jit^2*diag(nt))%*%z - yt)
    lam_collect = c(lam_collect,norm(lam,type='2'))
    
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    toc()
  }
  ploty(mc)
  l_fold1 = l
  p_fold1 = p
  
  #2 -fold redo the above
  xt = data$xv; xv = data$xt;yt = data$yv; yv = data$yt 
  l_collect = c();p_collect=c()
  
  lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')
  l = l_init = 0.9 
  p = p_init = 1.5
  l_collect = l
  p_collect = p
  l1 = abs(round(rnorm(1,mean = l_init,sd = 0.1),3))
  p1 = abs(round(rnorm(1,mean = p_init,sd = 0.1),3))
  
  z = z_init = solve(K(xt,xt,l1,p1)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
  
  gap_collect = norm(z-solve(K(xt,xt,l,p)+jit^2*diag(nt),yt),type="2")
  
  f_collect = c();p1_collect = c()
  
  cat("MC:",mc,"| l init:",green(l,p),"| jitter",green(l1,p1),
      blue(round(f(l,p,z,lam),4)),blue(round(partial1(l,p,z,lam),4)),'\n')
  
  for (i in 1:Maxiter){
    tic()
    ##==l update==##
    iter = 0
    while(norm(grad_l(l,p,z,lam),type="2")>1e-03){
      d = grad_l(l,p,z,lam)
      t = 0.125
      while(f(l-t*d,p,z,lam)>=f(l,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      l = l-t*d
      iter = iter+1
    }
    l_collect = c(l_collect,l)
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat(i,'l:',red(round(l,4)),yellow(iter),blue(round(f(l,p,z,lam),3)),'|')
    
    ##==p update==##
    iter = 0
    while(norm(grad_p(l,p,z,lam),type="2")>1e-03){
      d = grad_p(l,p,z,lam)
      t = 0.125
      while(f(l,p-t*d,z,lam)>=f(l,p,z,lam)-alpha*t*d*d){
        t = t*beta
      }
      p = p-t*d
      iter = iter+1
    }
    p_collect = c(p_collect,p)
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat('p:',red(round(p,4)),yellow(iter),blue(round(f(l,p,z,lam),3)),'|')
    
    ##==Z update==##
    iter = length(z)
    K_tt = K(xt,xt,l,p)
    K_vt = K(xv,xt,l,p)
    bv <- -2*crossprod(K_vt,yv) + 
      (K_tt+jit^2*diag(length(xt))) %*% lam - 
      rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
    sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
    
    d = -grad_z(l,p,z,lam)
    for (ii in 1:(iter-1)){
      stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
      zold = z
      z = z+stepsize*d
      if(norm(z-zold,type="2") <= 1e-5){break}
      bb = as.numeric((crossprod(grad_z(l,p,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
      d = -grad_z(l,p,z,lam)+bb*d
    }
    
    z_collect = c(z_collect,norm(z,type='2'))
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    cat('z norm:',red(norm(z,type='2')),yellow(ii),blue(round(f(l,p,z,lam),5)))
    
    gap = norm(z-solve(K(xt,xt,l,p)+jit^2*diag(nt),yt),type="2");gap_collect = c(gap_collect,gap)
    cat('| GAP',gap,'|',blue(round(f(l,p,z,lam),4)),'|')
    
    #dual update
    lam = lam + rhorho*((K(xt,xt,l,p)+jit^2*diag(nt))%*%z - yt)
    lam_collect = c(lam_collect,norm(lam,type='2'))
    
    f_collect = c(f_collect,f(l,p,z,lam))
    p1_collect = c(p1_collect,partial1(l,p,z,lam))
    toc()
  }
  ploty(mc)
  l_fold2 = l
  p_fold2 = p
  
  l_corpora = c(l_corpora,(l_fold1+l_fold2)/2)
  p_corpora = c(p_corpora,(p_fold1+p_fold2)/2)
}
#file = paste(format(Sys.time(), "%b_%d %X"),mc)
#write.csv(data.frame(l=l_corpora,p=p_corpora),file=paste(file,'_params.csv',sep=''),row.names = F)