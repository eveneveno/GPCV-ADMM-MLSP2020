#function collections
myfun_per <- function(x,y) {return(pi*abs(x-y))}
myfun_se <- function(x,y) {return((x-y)^2)}

K<- function(x,y,l1,l2,l3){
  inn_per = outer(x,y,myfun_per)
  inn_se = outer(x,y,myfun_se)
  return(theta1*exp(-inn_se/(2*l1^2)) + theta3*(exp(-inn_se/(2*l2^2))*exp(-2*(sin(inn_per/(0.197))^2)/(l3^2))))
}

Kg_l1 <- function(x,y,l1,l2,l3){
  inn_se = outer(x,y,myfun_se)
  return(theta1*exp(-inn_se/(2*l1^2)) * (inn_se/l1^3))
}

Kg_l2 <- function(x,y,l1,l2,l3){
  inn_se = outer(x,y,myfun_se)
  inn_per = outer(x,y,myfun_per)
  return(theta3*exp(-inn_se/(2*l2^2)) * (inn_se/l2^3) * exp(-2*(sin(inn_per/(0.197))^2)/(l3^2)))
}

Kg_l3 <- function(x,y,l1,l2,l3){
  inn_per = outer(x,y,myfun_per)
  inn_se = outer(x,y,myfun_se)
  return(theta3*exp(-2*(sin(inn_per/(0.197))^2)/(l3^2))* (4*(sin(inn_per/(0.197))^2)/(l3^3))*exp(-inn_se/(2*l2^2)))
}

f <- function(l1,l2,l3,z,lam){
  Cz = (K(xt,xt,l1,l2,l3)+jit^2*diag(nt))%*%z
  return(norm(yv-K(xv,xt,l1,l2,l3)%*%z,type="2")^2 + t(lam)%*%(Cz-yt) + rho/2*norm(Cz-yt,type="2")^2)
}

partial1 <- function(l1,l2,l3,z,lam){return(norm(yv-K(xv,xt,l1,l2,l3)%*%z,type="2")^2)}

grad_l1 <- function(l1,l2,l3,z,lam){
  Kg_tt = Kg_l1(xt,xt,l1,l2,l3); Kg_vt = Kg_l1(xv,xt,l1,l2,l3)
  K_tt = K(xt,xt,l1,l2,l3); K_vt = K(xv,xt,l1,l2,l3)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z +crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_l2 <- function(l1,l2,l3,z,lam){
  Kg_tt = Kg_l2(xt,xt,l1,l2,l3); Kg_vt = Kg_l2(xv,xt,l1,l2,l3)
  K_tt = K(xt,xt,l1,l2,l3); K_vt = K(xv,xt,l1,l2,l3)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z +crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_l3 <- function(l1,l2,l3,z,lam){
  Kg_tt = Kg_l3(xt,xt,l1,l2,l3); Kg_vt = Kg_l3(xv,xt,l1,l2,l3)
  K_tt = K(xt,xt,l1,l2,l3); K_vt = K(xv,xt,l1,l2,l3)
  
  grad = -2*crossprod(yv,Kg_vt)%*%z +crossprod(z,t(Kg_vt))%*%K_vt%*%z + crossprod(z,t(K_vt))%*%Kg_vt%*%z + 
    crossprod(lam,Kg_tt)%*%z + rho*crossprod(jit^2*z-yt,Kg_tt)%*%z + (rho/2)*(crossprod(z,K_tt)%*%Kg_tt%*%z + crossprod(z,Kg_tt)%*%K_tt%*%z)
  return(as.numeric(grad))
}

grad_z <- function(l1,l2,l3,z,lam){
  K_tt = K(xt,xt,l1,l2,l3)
  K_vt = K(xv,xt,l1,l2,l3)
  
  b <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt
  sigma <- 2*(crossprod(K_vt,K_vt%*%z) + 
                (rho/2)*(K_tt+jit^2*diag(length(xt))) %*% ((K_tt+jit^2*diag(length(xt)))%*%z))
  return(sigma + b)
}

ploty <- function(){
  par(mar=c(2,3,1,1))
  layout(matrix(c(1,2,3,4,5,6,7,7,7,7,7,7),4,3,byrow = TRUE))
  plot(x=0:(length(l1_collect)-1),y=l1_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l1')
  abline(h=l1_true,col=3,lty=2)
  abline(h=l1,col=4)
  text(round(l1,5),x=(length(l1_collect)-20), y=l1+0.015)
  
  plot(x=0:(length(l2_collect)-1),y=l2_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l2')
  abline(h=l2_true,col=3,lty=2)
  abline(h=l2,col=4)
  text(round(l2,5),x=(length(l2_collect)-20), y=l2+0.015)
  
  plot(x=0:(length(l3_collect)-1),y=l3_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of l3')
  abline(h=l3_true,col=3,lty=2)
  abline(h=l3,col=4)
  text(round(l3,5),x=(length(l3_collect)-20), y=l3+0.015)
  
  plot(x=0:(length(z_collect)-1),y=z_collect,type='l',col=2,
       xlab = 'Iter',ylab='change of z (l2 norm)')
  abline(h=norm(z,type="2"),col=4)
  text(round(norm(z,type="2"),6),x=(length(z_collect)-20), y=norm(z,type="2")+0.015) 
  
  plot(x=0:(length(gap_collect)-1),y=gap_collect,type='l',col=2,
       xlab = 'Iter',ylab='Inequality Gap Control')
  text(paste('Rho:',rho),x=(length(gap_collect)/2), y=(max(gap_collect)+min(gap_collect))/1.8) 
  text(paste('RhoRho:',rhorho),x=(length(gap_collect)/2), y=(max(gap_collect)+min(gap_collect))/2.2) 
  
  plot(p1_collect,type='l',cex=0.5,col=2,ylab='Original objective')
  abline(h=tail(p1_collect, n=1),col=4) 
  text(round(tail(p1_collect,n=1),6),x=(length(p1_collect)-20), y=tail(p1_collect, n=1)+0.015) 
  
  plot(xt,yt,type='p',col='darkgrey',xlim=c(0,12.5),ylim=c(-2,3),
       cex=1,xlab='x',ylab='CO2 concentration',xaxt='n',yaxt='n',pch=16)
  ticks = x[seq(25,length(x),120)]
  axis(1,at= c(ticks,tail(ticks,n=1)+1.97,tail(ticks,n=1)+2*1.97),labels=c(1960,1970,1980,1990,2000,2010,2020),font=2)
  points(xv,yv,cex=1,pch=16)
}

ploty2 <- function(){
  z_true = solve(K(x,x,l1_true,l2_true,l3_true)+jit^2*diag(length(x)),y)
  y_pred = K(xv,x,l1_true,l2_true,l3_true)%*%z_true #(2) TRUE
  cat(norm(yv-y_pred,type="2"))
  lines(xv,y_pred,col='green',lwd=2)
  y_pred_date = K(pred_date,x,l1_true,l2_true,l3_true)%*%z_true #variance not shown
  lines(pred_date,y_pred_date,col='green',type='l',lty=2,lwd=2)
  
  z_true = solve(K(x,x,l1,l2,l3)+jit^2*diag(length(x)),y)
  y_pred = K(xv,x,l1,l2,l3)%*%z_true #(2) TRUE
  cat(norm(yv-y_pred,type="2"))
  lines(xv,y_pred,col='red',lwd=2)
  y_pred_date = K(pred_date,x,l1,l2,l3)%*%z_true #variance not shown
  lines(pred_date,y_pred_date,col='red',type='l',lty=2,lwd=2)
  
  z_true = solve(K(x,x,l1_init,l2_init,l3_init)+jit^2*diag(length(x)),y)
  y_pred = K(xv,x,l1_init,l2_init,l3_init)%*%z_true #(2) TRUE
  cat(norm(yv-y_pred,type="2"))
  lines(xv,y_pred,col='blue',lwd=2)
  y_pred_date = K(pred_date,x,l1_init,l2_init,l3_init)%*%z_true #variance not shown
  lines(pred_date,y_pred_date,col='blue',type='l',lty=2,lwd=2)
}

#parameter setting
jit = 0.5
alpha = beta = 0.5
rho = 5
rhorho = 5
Maxiter = 50

#rescaled output variance
theta1 = 66*0.01316*4
theta3 = 2.4*0.01316*4

library(MASS)
library(tictoc)
library(crayon)