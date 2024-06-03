library(INLA)
inla.setOption(inla.mode="experimental")
inla.setOption(num.threads=10) 

n2<-200 # number of subjects
n3<-20 # number of families
m<-10 # maximum number of observed times for each subjects    

n.iteration<-500

for (l in 1:n.iteration){
  
  # importing data
  data<-read.csv(paste0("data_",l,".txt"))
  
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
  
  l.long <- c(z, rep(NA, n4), rep(NA, n2))
  z.long <- c(rep(NA, n1), data.l$y, rep(NA, n2))
  y.term <- inla.surv(time = c(rep(NA, n1+n4), time.D), event = c(rep(NA, n1+n4), delta.D.fi))
  y.joint <- list(l.long, z.long, y.term)
  
  linear.covariate <- data.frame(mu = c(rep(NA,n1+n4),rep(1,n2)),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 l.T.R = c(log(data$T.R), rep(0, n4), rep(0, n2)),  # related to logit(pi): logistic component for longitudinal outcomes
                                 intercept.l = c(rep(1,n1), rep(0, n4), rep(0,n2)),
                                 x1.l = c(data$x1.l, rep(0, n4), rep(0, n2)), 
                                 x2.l = c(data$x2.l, rep(0, n4), rep(0, n2)),
                                 z.T.R = c(rep(0, n1), log(data.l$T.R), rep(0, n2)),  # related to log(mu): negative binomial component for longitudinal outcomes
                                 intercept.z = c(rep(0,n1), rep(1,n4), rep(0,n2)),
                                 x1.z = c(rep(0, n1), data.l$x1.l, rep(0, n2)),
                                 x2.z = c(rep(0, n1), data.l$x2.l, rep(0, n2)),
                                 x1.d = c(rep(0, n1+n4), data.d$x1.d), # related to hazard function for terminal event
                                 x2.d = c(rep(0, n1+n4), data.d$x2.d))
  
  random.covariate <- list(l.YL = c(data$id, rep(NA, n4), rep(NA, n2)), 
                           l.ZL = c(data$famid, rep(NA, n4), rep(NA, n2)), 
                           z.YL = c(rep(NA, n1), data.l$id, rep(NA, n2)), 
                           z.ZL = c(rep(NA, n1), data.l$famid, rep(NA, n2)), 
                           d.YL = c(rep(NA, n1+n4), (1:n2)), 
                           d.ZL = c(rep(NA, n1+n4), data.d$famid))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  formula = Y ~ - 1 + mu + intercept.l + x1.l + x2.l + offset(l.T.R) +
    intercept.z + x1.z + x2.z + offset(z.T.R) +
    x1.d + x2.d + 
    f(d.YL, model="iid", n=n2, hyper=list(prec = list(param = c(2.5, 0.5)))) +
    f(l.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(z.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(d.ZL, model="iid", n=n3, hyper=list(prec = list(param = c(2.5, 0.5)))) +
    f(l.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(z.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) 
  
  inla.model <- inla(formula, family = c("binomial","zeroinflatednbinomial0","weibullsurv"),
                     data = data.f,control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE,openmp.strategy="huge"),
                     control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE),size = list(param=3))),list(variant=1)),
                     control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE, control.fixed = list(prec.intercept = 0.1))
  
  inla.model.result<-summary(inla.model)
  result<-capture.output(inla.model.result)
  file = file(paste0("result_",l,".txt"))
  writeLines(result, file)
  close(file)
}
