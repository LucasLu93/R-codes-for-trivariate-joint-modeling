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
  
  y.recu <- inla.surv(time = c(time.R, rep(NA, n2)), event = c(delta.R.fij, rep(NA, n2)))
  y.term <- inla.surv(time = c(rep(NA, n1), time.D), event = c(rep(NA, n1), delta.D.fi))
  y.joint <- list(y.recu, y.term)
  
  linear.covariate <- data.frame(mu = as.factor(c(rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                                 x1.r = c(data$x1.r, rep(0, n2)), # relative to hazard function for recurrent events
                                 x2.r = c(data$x2.r, rep(0, n2)),
                                 x1.d = c(rep(0, n1), data.d$x1.d), # related to hazard function for terminal event
                                 x2.d = c(rep(0, n1), data.d$x2.d))
  
  random.covariate <- list(r.YR = c(data$id, rep(NA, n2)), 
                           r.ZR = c(data$famid, rep(NA, n2)), 
                           d.YR = c(rep(NA, n1), 1:n2), 
                           d.ZR = c(rep(NA, n1), data.d$famid))
  
  data.f <- c(linear.covariate,random.covariate)
  data.f$Y <- y.joint
  
  formula = Y ~ - 1 + mu + 
    x1.r + x2.r + 
    x1.d + x2.d + 
    f(d.YR, model="iid", n=n2, hyper=list(prec = list(param = c(2.5, 0.5)))) +
    f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
    f(r.ZR, model="iid", n=n3, hyper=list(prec = list(param = c(2.5, 0.5)))) +
    f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

  inla.model <- inla(formula, family = c("weibullsurv","weibullsurv"),
                     data = data.f,control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE,openmp.strategy="huge"),
                     control.family = list(list(variant=1),list(variant=1)),
                     control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE, control.fixed = list(prec.intercept = 0.1))
  
  inla.model.result<-summary(inla.model)
  result<-capture.output(inla.model.result)
  file = file(paste0("result_",l,".txt"))
  writeLines(result, file)
  close(file)
}
  
  
  
  