source('preparation.r')
source('read.r')
options(digits = 4)

# Rescale the data
x = 10*(x-min(x))/(max(x)-min(x))
y = 4*(y-min(y))/(max(y)-min(y)) 
y = y-mean(y)

# Data plot
# plot(x,y,xlim=c(0,15),ylim=c(-2,4),type='l')

# Random split data to 2-folds
train_idx = sort(sample(1:length(x),length(x)*1/2))

xt = x[train_idx];yt = y[train_idx];nt = length(xt)
xv = x[-train_idx];yv = y[-train_idx]
# points(xt,yt,col='gray',cex=0.2)
# points(xv,yv,col='black',cex=0.2)

# Data for prediction
pred_date = 10+x[1:140]

# Rescaled estimates given by [Rasmussen and Williams, 2006]
# l1_true = 13.2 #67*0.0197*10
# l2_true = 17.73 #90*0.0197*10
# l3_true = 0.2955 #1.5*0.0197*10

# Initialization of Hyper-parameter 
# l1 = l1_init = 30
# l2 = l2_init = 10
# l3 = l3_init = 0.5

# l1 = l1_init = 100
# l2 = l2_init = 15
# l3 = l3_init = 1

# Perturbed initialization of Auxiliary variable
l11 = abs(round(rnorm(1,mean = l1_init,sd = 0.1),3))
l21 = abs(round(rnorm(1,mean = l2_init,sd = 0.1),3))
l31 = abs(round(rnorm(1,mean = l3_init,sd = 0.1),3))
z = z_init = solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')

# Collect the change of variables
l1_collect = l1; l2_collect = l2; l3_collect = l3
cost_collect = c();part1_collect = c() 
gap_collect = norm(z-solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt),type="2")
f_collect = c();p1_collect = c()

# ADMM updating scheme
# K = 1 
for (i in 1:Maxiter){
  if (i == 1){
    cat("l1 l2 l3 init:",green(l1,l2,l3),"| jitter",green(l11,l21,l31),
        blue(round(f(l1,l2,l3,z,lam),4)),blue(round(partial1(l1,l2,l3,z,lam),4)),'\n') 
  }
  tic()
  ##==l1 update==##
  iter = 0
  while(norm(grad_l1(l1,l2,l3,z,lam),type="2")>1e-03){
    d = grad_l1(l1,l2,l3,z,lam)
    t = 1
    while(f(l1-t*d,l2,l3,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l1 = l1-t*d
    iter = iter+1
  }
  l1_collect = c(l1_collect,l1)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l1,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')
  
  ##==l2 update==##
  iter = 0
  while(norm(grad_l2(l1,l2,l3,z,lam),type="2")>1e-03){
    d = grad_l2(l1,l2,l3,z,lam)
    t = 1
    while(f(l1,l2-t*d,l3,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l2 = l2-t*d
    iter = iter+1
  }
  l2_collect = c(l2_collect,l2)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l2,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')  
  
  ##==l3 update==##
  iter = 0
  while(norm(grad_l3(l1,l2,l3,z,lam),type="2")>1e-02){
    d = grad_l3(l1,l2,l3,z,lam)
    t = 1
    while(f(l1,l2,l3-t*d,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l3 = l3-t*d
    iter = iter+1
  }
  l3_collect = c(l3_collect,l3)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l3,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')  
  
  ##==Z update==##
  iter = length(z)
  K_tt = K(xt,xt,l1,l2,l3)
  K_vt = K(xv,xt,l1,l2,l3)
  bv <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
  sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
  
  d = -grad_z(l1,l2,l3,z,lam)
  for (ii in 1:(iter-1)){
    stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
    zold = z
    z = z+stepsize*d
    if(norm(z-zold,type="2") <= 1e-5){break}
    bb = as.numeric((crossprod(grad_z(l1,l2,l3,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
    d = -grad_z(l1,l2,l3,z,lam)+bb*d
  }
  z_collect = c(z_collect,norm(z,type='2'))
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(norm(z,type='2'),yellow(ii),blue(round(f(l1,l2,l3,z,lam),5)))
  
  # compute the inequality gap
  gap = norm(z-solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt),type="2")
  gap_collect = c(gap_collect,gap)
  cat('gap',gap,' |',blue(round(f(l1,l2,l3,z,lam),4)),' |')
  
  lam = lam + rhorho*((K(xt,xt,l1,l2,l3)+jit^2*diag(nt))%*%z - yt)
  lam_collect = c(lam_collect,norm(lam,type='2'))
  cat(blue(round(f(l1,l2,l3,z,lam),4)),' |')
  
  toc()
}

l1_1 = l1; l2_1 = l2; l3_1 = l3

# K = 2 (exchange training & validation data set)
# redo the above procedure
# Initialization of Hyper-parameter 
# l1 = l1_init = 30
# l2 = l2_init = 10
# l3 = l3_init = 0.5

# l1 = l1_init = 100
# l2 = l2_init = 15
# l3 = l3_init = 1

# Perturbed initialization of Auxiliary variable
l11 = abs(round(rnorm(1,mean = l1_init,sd = 0.1),3))
l21 = abs(round(rnorm(1,mean = l2_init,sd = 0.1),3))
l31 = abs(round(rnorm(1,mean = l3_init,sd = 0.1),3))
z = z_init = solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt);z_collect = norm(z,type='2')
lam = matrix(1,nt,1);lam_collect = norm(lam,type='2')

# Collect the change of variables
l1_collect = l1; l2_collect = l2; l3_collect = l3
cost_collect = c();part1_collect = c() 
gap_collect = norm(z-solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt),type="2")
f_collect = c();p1_collect = c()

xv = x[train_idx];yv = y[train_idx];nt = length(xt)
xt = x[-train_idx];yt = y[-train_idx]

for (i in 1:Maxiter){
  if (i == 1){
    cat("l1 l2 l3 init:",green(l1,l2,l3),"| jitter",green(l11,l21,l31),
        blue(round(f(l1,l2,l3,z,lam),4)),blue(round(partial1(l1,l2,l3,z,lam),4)),'\n') 
  }
  tic()
  ##==l1 update==##
  iter = 0
  while(norm(grad_l1(l1,l2,l3,z,lam),type="2")>1e-03){
    d = grad_l1(l1,l2,l3,z,lam)
    t = 1
    while(f(l1-t*d,l2,l3,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l1 = l1-t*d
    iter = iter+1
  }
  l1_collect = c(l1_collect,l1)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l1,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')
  
  ##==l2 update==##
  iter = 0
  while(norm(grad_l2(l1,l2,l3,z,lam),type="2")>1e-03){
    d = grad_l2(l1,l2,l3,z,lam)
    t = 1
    while(f(l1,l2-t*d,l3,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l2 = l2-t*d
    iter = iter+1
  }
  l2_collect = c(l2_collect,l2)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l2,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')  
  
  ##==l3 update==##
  iter = 0
  while(norm(grad_l3(l1,l2,l3,z,lam),type="2")>1e-02){
    d = grad_l3(l1,l2,l3,z,lam)
    t = 1
    while(f(l1,l2,l3-t*d,z,lam)>=f(l1,l2,l3,z,lam)-alpha*t*d*d){
      t = t*beta
    }
    l3 = l3-t*d
    iter = iter+1
  }
  l3_collect = c(l3_collect,l3)
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(i,red(round(l3,4)),yellow(iter),blue(round(f(l1,l2,l3,z,lam),3)),'|')  
  
  ##==Z update==##
  iter = length(z)
  K_tt = K(xt,xt,l1,l2,l3)
  K_vt = K(xv,xt,l1,l2,l3)
  bv <- -2*crossprod(K_vt,yv) + 
    (K_tt+jit^2*diag(length(xt))) %*% lam - 
    rho* (K_tt+jit^2*diag(length(xt))) %*% yt  
  sigma <- 2*(crossprod(K_vt) + (rho/2)*crossprod(K_tt+jit^2*diag(length(xt))))
  
  d = -grad_z(l1,l2,l3,z,lam)
  for (ii in 1:(iter-1)){
    stepsize = -as.numeric(crossprod(d,sigma%*%z+bv)/(crossprod(d,sigma)%*%d))
    zold = z
    z = z+stepsize*d
    if(norm(z-zold,type="2") <= 1e-5){break}
    bb = as.numeric((crossprod(grad_z(l1,l2,l3,z,lam),sigma)%*%d)/crossprod(d,sigma)%*%d)
    d = -grad_z(l1,l2,l3,z,lam)+bb*d
  }
  z_collect = c(z_collect,norm(z,type='2'))
  f_collect = c(f_collect,f(l1,l2,l3,z,lam))
  p1_collect = c(p1_collect,partial1(l1,l2,l3,z,lam))
  cat(norm(z,type='2'),yellow(ii),blue(round(f(l1,l2,l3,z,lam),5)))
  
  # compute the inequality gap
  gap = norm(z-solve(K(xt,xt,l1,l2,l3)+jit^2*diag(nt),yt),type="2")
  gap_collect = c(gap_collect,gap)
  cat('gap',gap,' |',blue(round(f(l1,l2,l3,z,lam),4)),' |')
  
  lam = lam + rhorho*((K(xt,xt,l1,l2,l3)+jit^2*diag(nt))%*%z - yt)
  lam_collect = c(lam_collect,norm(lam,type='2'))
  cat(blue(round(f(l1,l2,l3,z,lam),4)),' |')
  
  toc()
}

l1_2 = l1; l2_2 = l2; l3_2 = l3

l1 = (l1_1+l1_2)/2
l2 = (l2_1+l2_2)/2
l3 = (l3_1+l3_2)/2

# ploty()
# ploty2()
