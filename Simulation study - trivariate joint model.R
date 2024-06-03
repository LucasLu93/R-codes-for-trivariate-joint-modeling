library(mvtnorm)
library(distributions3)
library(countreg)
library(INLA)
library(IDPmisc)
inla.setOption(inla.mode="experimental")
inla.setOption(num.threads=10) 

n2<-200 # number of subjects
n3<-20 # number of families
m<-10 # maximum number of observed times for each subjects    
id<-1:n2
id.2<-rep(id,each=m)

# random effects
sd.Y.R<-1.5 # related to recurrent events
var.Y.R<-sd.Y.R^2
sd.Y.L<-1.8 # related to longitudinal outcomes
var.Y.L<-sd.Y.L^2
cor.Y<-0.9
cov.Y<-cor.Y*sqrt(var.Y.R*var.Y.L)
vcov.Y<-matrix(c(var.Y.R,cov.Y,cov.Y,var.Y.L),2,2)

sd.Z.R<-0.5 # related to recurrent events
var.Z.R<-sd.Z.R^2
sd.Z.L<-0.7 # related to longitudinal outcomes
var.Z.L<-sd.Z.L^2
cor.Z<-0.4
cov.Z<-cor.Z*sqrt(var.Z.R*var.Z.L)
vcov.Z<-matrix(c(var.Z.R,cov.Z,cov.Z,var.Z.L),2,2)

# terminal event
beta<-c(3.5,-1.1)
xi<--0.5
lambda.d<-0.7 # shape for weibull distribution
rho.d<-33.5 # scale for weibull distribution

# recurrent event
gamma<-c(0.3,-0.3)
alpha<-0.4
lambda.r<-1.2 # shape for weibull distribution
rho.r<-0.2 # scale for weibull distribution

# longitudinal outcomes
chi<-c(1.6,-0.3,0.02)
theta<-c(-1.9,1.6,-1.1)
omega.l<--0.6
omega.z<-0.7
tau.l<-1.2
tau.z<-1.2
varphi<-8.5 # dispersion parameter

n.iteration<-500

for (l in 1:n.iteration){
  
  famid<-rep(1:n3,each=10)
  famid.2<-rep(famid,each=m)
    
  # random effects
  Y<-rmvnorm(n2,c(0,0),vcov.Y)
  Z<-rmvnorm(n3,c(0,0),vcov.Z)
    
  Y.2<-Y[id.2,] # subject-specific random effect
  Z.2<-Z[famid.2,] # family-specific random effect
    
  # terminal event
  x1<-rbinom(n2,1,0.4) 
  x1.d<-x1
  x2<-rnorm(n3,-0.2,0.5) 
  x2.d<-x2[famid] 
    
  mu.d<-apply(cbind(x1.d,x2.d,Y[,2],Z[famid,][,2],Y[,1],Z[famid,][,1]),1,function(x){
    beta[1]*x[1]+beta[2]*x[2]+x[3]+x[4]+x[5]+xi*x[6]})
  surv.d<-sapply(mu.d,function(x){rweibull(1,lambda.d,rho.d/exp(x))}) 
    
  q<-runif(n2,0,0.6) 
  status.d<-rep(1,n2) 
  for (i in 1:n2){
    if (surv.d[i]>q[i]) {
      surv.d[i]<-q[i]
      status.d[i]<-0
    }
  }
    
  x1.d<-rep(x1.d,each=m)
  x2.d<-rep(x2.d,each=m)
  surv.d<-rep(surv.d,each=m)
  status.d<-rep(status.d,each=m)
    
  # recurrent event
  x1.r<-x1[id.2]
  x2.r<-x2[famid.2] 
    
  mu.r<-apply(cbind(x1.r,x2.r,Y.2[,1],Z.2[,1]),1,function(x){
    gamma[1]*x[1]+gamma[2]*x[2]+alpha*x[3]+x[4]})
  surv.r<-sapply(mu.r,function(x){rweibull(1,lambda.r,rho.r/exp(x))}) # also the gap time between two consecutive longitudinal outcomes
    
  cum.gap<-surv.r[1] # cumulative gap time for recurrent events
  for (i in 1:n2){
    for (j in 1:m){
      if (j==1) (cum.gap[(i-1)*m+j]<-surv.r[(i-1)*m+j])
      if (j>1) (cum.gap[(i-1)*m+j]<-cum.gap[(i-1)*m+j-1]+surv.r[(i-1)*m+j])
    }
  }
    
  status.r<-rep(1,n2*m) 
  v<-rep(1,n2*m) # indicator variable for valid data
  for (i in 1:n2){
    for (j in 1:m){
      if (j==1 & cum.gap[(i-1)*m+1]>surv.d[(i-1)*m+1]) {
        v[((i-1)*m+2):((i-1)*m+m)]<-0
        surv.r[(i-1)*m+1]<-surv.d[(i-1)*m+1]
        status.r[(i-1)*m+1]<-0
      }
      if (j>1 & j<m & cum.gap[(i-1)*m+j]>surv.d[(i-1)*m+j]) {
        v[((i-1)*m+j+1):((i-1)*m+m)]<-0
        surv.r[(i-1)*m+j]<-surv.d[(i-1)*m+j]-cum.gap[(i-1)*m+j-1]
        status.r[(i-1)*m+j]<-0
      }
      if (j==m & cum.gap[(i-1)*m+m]>surv.d[(i-1)*m+m]) {
        surv.r[(i-1)*m+m]<-surv.d[(i-1)*m+m]-cum.gap[(i-1)*m+m-1]
        status.r[(i-1)*m+m]<-0
      }
    }
  }
    
  # longitudinal outcomes
  T.R<-surv.r 
  log.T.R<-log(T.R) # offset term
  log.T.R[is.na(log.T.R)]<-0
  x1.l<-x1[id.2]
  x2.l<-x2[famid.2] 
    
  # logistic component
  mu.l<-apply(cbind(log.T.R,x1.l,x2.l,Y.2[,2],Z.2[,2]),1,function(x){
    x[1]+chi[1]+chi[2]*x[2]+chi[3]*x[3]+omega.l*x[4]+tau.l*x[5]}) 
  z<-sapply(mu.l,function(x){rbinom(1,1,exp(x)/(1+exp(x)))})
  
  # negative binomial component
  mu.z<-apply(cbind(log.T.R,x1.l,x2.l,Y.2[,2],Z.2[,2]),1,function(x){
    x[1]+theta[1]+theta[2]*x[2]+theta[3]*x[3]+omega.z*x[4]+tau.z*x[5]}) 
    
  y<-NULL # longitudinal outcomes
  for (i in 1:(n2*m)){
    if (z[i]==1) (y[i]<-0)
    if (z[i]==0) (y[i]<-rztnbinom(1,mu=exp(mu.z[i]),theta=varphi))
  }
    
  # generating data
  data<-cbind(id.2,famid.2,T.R,x1.l,x2.l,x1.r,x2.r,surv.r,status.r,x1.d,x2.d,surv.d,status.d,y)
  data<-data[v==1,]
  d<-NULL # indicator variable of each subject for terminal events
  n<-nrow(data)
  for (i in 1:(n-1)){
    if (data[i,1]==data[i+1,1]) {
      d[i]<-0} else {
        d[i]<-1}}
  d[n]<-1
  data<-data.frame(cbind(data,d))
  write.table(data,paste0("data_",l,".txt"))
    
  ### Constructing Model
  # longitudinal outcomes
  data.l<-data[data$y!=0,]
  n4<-nrow(data.l)
    
  # preparation for recurrent events
  time.R<-data$surv.r 
  delta.R.fij<-data$status.r 
    
  # preparation for terminal event
  data.d<-data[data$d==1,]
  time.D<-data.d$surv.d 
  delta.D.fi<-data.d$status.d 
    
  # other preparations
  n1<-nrow(data) # number of observations
  z<-rep(0,n1)
  z[data$y==0]<-1
    
  l.long <- c(z, rep(NA, n4),  rep(NA, n1), rep(NA, n2))
  z.long <- c(rep(NA, n1), data.l$y, rep(NA, n1), rep(NA, n2))
  y.recu <- inla.surv(time = c(rep(NA, n1+n4), time.R, rep(NA, n2)), event = c(rep(NA,n1+n4), delta.R.fij, rep(NA, n2)))
  y.term <- inla.surv(time = c(rep(NA, n1+n4+n1), time.D), event = c(rep(NA, n1+n4+n1), delta.D.fi))
  y.joint <- list(l.long, z.long, y.recu, y.term)
    
  linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n1+n4),rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.T.R = c(log(data$T.R), rep(0, n1+n4), rep(0, n2)),  # related to logit(pi): logistic component for longitudinal outcomes
                                 intercept.l = c(rep(1,n1), rep(0, n1+n4), rep(0,n2)),
                                 x1.l = c(data$x1.l, rep(0, n1+n4), rep(0, n2)), 
                                 x2.l = c(data$x2.l, rep(0, n1+n4), rep(0, n2)),
                                 z.T.R = c(rep(0, n1), log(data.l$T.R), rep(0,n1), rep(0, n2)),  # related to log(mu): negative binomial component for longitudinal outcomes
                                 intercept.z = c(rep(0,n1), rep(1,n4), rep(0,n1), rep(0,n2)),
                                 x1.z = c(rep(0, n1), data.l$x1.l, rep(0,n1), rep(0, n2)),
                                 x2.z = c(rep(0, n1), data.l$x2.l, rep(0,n1), rep(0, n2)),
                                 x1.r = c(rep(0, n1+n4), data$x1.r, rep(0, n2)), # relative to hazard function for recurrent events
                                 x2.r = c(rep(0, n1+n4), data$x2.r, rep(0, n2)),
                                 x1.d = c(rep(0, n1+n4+n1), data.d$x1.d), # related to hazard function for terminal event
                                 x2.d = c(rep(0, n1+n4+n1), data.d$x2.d))
    
  random.covariate <- list(l.YL = c(data$id+n2, rep(NA, n1+n4), rep(NA, n2)), 
                           l.ZL = c(data$famid+n3, rep(NA, n1+n4), rep(NA, n2)), 
                           z.YL = c(rep(NA, n1), data.l$id+n2, rep(NA, n1), rep(NA, n2)), 
                           z.ZL = c(rep(NA, n1), data.l$famid+n3, rep(NA, n1), rep(NA, n2)), 
                           r.YR = c(rep(NA, n1+n4), data$id, rep(NA, n2)), 
                           r.ZR = c(rep(NA, n1+n4), data$famid, rep(NA, n2)), 
                           d.YR = c(rep(NA, n1+n4+n1), 1:n2), 
                           d.ZR = c(rep(NA, n1+n4+n1), data.d$famid),
                           d.YL = c(rep(NA, n1+n4+n1), (1:n2)+n2), 
                           d.ZL = c(rep(NA, n1+n4+n1), data.d$famid+n3))
    
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
    
  formula = Y ~ - 1 + mu + intercept.l + x1.l + x2.l + offset(l.T.R) +
    intercept.z + x1.z + x2.z + offset(z.T.R) +
    x1.r + x2.r + 
    x1.d + x2.d + 
    f(d.YR, model="iidkd", order=2, n=2*n2, hyper = list(theta1 = list(param = c(5, 1, 1, 0)))) +
    f(d.YL, copy="d.YR") +
    f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(l.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(z.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(r.ZR, model="iidkd", order=2, n=2*n3, hyper = list(theta1 = list(param = c(5, 1, 1, 0)))) +
    f(d.ZL, copy="r.ZR") +
    f(l.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(z.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

  inla.model<-try(inla(formula, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                       data = data.f,control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE,openmp.strategy="huge"),
                       control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE),size = list(param=3))),list(variant=1),list(variant=1)),
                       control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE, control.fixed = list(prec.intercept = 0.1)), 
                  silent=TRUE)
  
  inla.model.result<-summary(inla.model)
  result<-capture.output(inla.model.result)
  file = file(paste0("result_",l,".txt"))
  writeLines(result, file)
  close(file)
}




