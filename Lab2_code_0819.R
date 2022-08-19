#Lab 2 code
##1
library('LaplacesDemon')
library('MASS')
library('mvtnorm')
data<-read.table(file='TempLambohov.txt',header = TRUE)
time<-as.vector(data$time)
temp<-as.vector(data$temp)
mu0<-c(-10,100,-100)
omega0<-0.02*diag(3)
v0<-3
sigma2<-2
#simulating
ndraws<-20
beta_simul<-matrix(1,ndraws,3)
for(i in 1:ndraws){
  sigma2_simul<-rinvchisq(1,v0,sigma2)
  beta_simul[i,]<-mvrnorm(1,mu0,sigma2_simul*solve(omega0))
}
plot(x=time,y=temp,xlim = c(0,1),ylim = c(-50,50),type = 'p',col='blue',main = 'mu=(-10,100,-100)')
for(i in 1:ndraws){
  b0<-beta_simul[i,1]
  b1<-beta_simul[i,2]
  b2<-beta_simul[i,3]
  curve(b0+b1*x+b2*x**2,0,1,col='red',add=TRUE)
}
#simulating from posterior
X<-as.matrix(data.frame(cons=rep(1,length=365),x=time,x2=time**2))
beta_h<-solve(t(X)%*%X)%*%t(X)%*%temp
mu_n<-solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%beta_h+omega0%*%mu0)
omega_n<-t(X)%*%X+omega0
v_n<-v0+nrow(data)
sigma2_n<-(v0*sigma2+(t(temp)%*%temp+t(mu0)%*%omega0%*%mu0-t(mu_n)%*%omega_n%*%mu_n))/v_n
sigma_pos<-vector(length = 1000)
beta_pos<-matrix(1,1000,3)
for(i in 1:1000){
  sigma_pos[i]<-rinvchisq(1,v_n,sigma2_n)
  beta_pos[i,]<-mvrnorm(1,mu_n,sigma_pos[i]**2*solve(omega_n))
}
hist(sigma_pos,main = 'histogram for sigma2',xlab = 'sigma2',breaks = 50)
hist(beta_pos[,1],main = 'histogram for beta0',xlab = 'beta0',breaks = 50)
hist(beta_pos[,2],main = 'histogram for beta1',xlab = 'beta1',breaks = 50)
hist(beta_pos[,3],main = 'histogram for beta2',xlab = 'beta2',breaks = 50)
plot(time,temp,type = 'p',col='blue',ylim = c(-30,30))
lines(time,apply(X%*%t(beta_pos),1,median),col='green')
lines(time,apply(X%*%t(beta_pos),1,quantile,prob=0.025),col='red')
lines(time,apply(X%*%t(beta_pos),1,quantile,prob=0.975),col='red')
x_tilde<--beta_pos[,2]/(2*beta_pos[,3])
hist(x_tilde,main = 'Density of x_tilde',breaks = 50)
data<-read.table("WomenAtWork.dat",header=TRUE)
X<-as.matrix(data[,2:8])
y<-as.matrix(data[,1])
##2
mcov<-diag(25,7)
logpost<-function(b){
  loglike<-sum((X%*%b)*y-log(1+exp(X%*%b)))
  if (abs(loglike) == Inf){logLik = -10000}
  logPrior<-dmvnorm(b,rep(0,7),mcov,log=TRUE)
  return(loglike+logPrior)
}
b_opt<-optim(rnorm(7,0,1),fn=logpost,method="BFGS",control=list(fnscale=-1),hessian=TRUE)
b_est<-b_opt$par
J<-solve(-b_opt$hessian)
c(b_est[6]-1.96*sqrt(J[6,6]),b_est[6]+1.96*sqrt(J[6,6]))
glmModel<-glm(Work ~ 0 + ., data = data, family = binomial)
summary(glmModel)
x_new<-c(1,20,12,8,43,7,10)
simul<-function(n){
  beta<-rmvnorm(n,mean=b_est,sigma=J)
  p<-apply(beta,1,FUN=function(b){1/(1+exp(x_new%*%b))})
  return(p)
}
h<-hist(simul(1000),breaks=10,plot=FALSE)
h$counts=h$counts/sum(h$counts)
plot(h)
x_new<-c(1,20,12,8,43,7,10)
simul_bi<-function(n){
  p<-simul(n)
  n_simul<-rbinom(1000,size=11,prob=p)
  return(n_simul)
}
h<-hist(simul_bi(1000),plot=FALSE)
h$counts=h$counts/sum(h$counts)
plot(h,main='posterior predictive distribution for the number of working women')
