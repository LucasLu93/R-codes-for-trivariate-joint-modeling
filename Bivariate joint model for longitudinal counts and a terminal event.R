rm(list=ls())
library(INLA)
inla.setOption(inla.mode="experimental")
inla.setOption(num.threads=1) 

# import the data
data<-read.table("~/Desktop/NFdat_newerdate.txt", stringsAsFactors=FALSE)

# preparation for longitudinal outcomes
# remove the observations with num >= 99
data$num[data$num>98] = NA

n_1<-nrow(data)-1
for (i in 1:n_1){
  if (!is.na(data$polyps[i]) & data$polyps[i] == 1 & is.na(data$num[i])) {
    data$num[i] = 1
  }
  if (!is.na(data$polyps[i]) & data$polyps[i] == 0 & is.na(data$num[i])) {
    data$num[i] = 0
  }
  if (is.na(data$polyps[i]) & is.na(data$num[i]) & data$id[i] == data$id[i+1] & data$polyps0[i+1] == 1) {
    data$num[i] = 1
    data$polyps[i] = 1
  }
  if (is.na(data$polyps[i]) & is.na(data$num[i]) & data$id[i] == data$id[i+1] & data$polyps0[i+1] == 0) {
    data$num[i] = 0
    data$polyps[i] = 0
  }
}

data.long<-data[!is.na(data$num),]
z<-rep(0,nrow(data.long)) # 282
z[data.long$num==0]<-1

data.long.z<-data.long[data.long$num!=0,]

# preparation for recurrent events
time.R<-data$gap
delta.R.fij<-data$visit

# preparation for terminal event 
data.term<-data[data$type=="crc",]
time.D<-data.term$crcgap 
delta.D.fi<-data.term$crc1 

# other preparations
n1<-nrow(data) # 459
n2<-nrow(data.term) # 177
n3<-length(unique(data.term$famid)) # 18
n4<-nrow(data.long) # 282
n5<-nrow(data.long.z) # 80

l.long <- c(z, rep(NA, n5), rep(NA, n2))
z.long <- c(rep(NA, n4), data.long.z$num, rep(NA, n2))
y.term <- inla.surv(time = c(rep(NA, n4+n5), time.D/max(time.D)), event = c(rep(NA, n4+n5), delta.D.fi))
y.joint <- list(l.long, z.long, y.term)

linear.covariate <- data.frame(mu = c(rep(NA,n4+n5),rep(1,n2)),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                               l.Time = c(log(data.long$gap), rep(0,n5+n2)), # related to logit(pi): logistic component for longitudinal outcomes 
                               intercept.l = c(rep(1,n4), rep(0,n5+n2)),
                               l.Sex = c(data.long$sex01, rep(0,n5+n2)), 
                               l.Age = c(data.long$age00_trans, rep(0,n5+n2)), 
                               l.Pcrc1age = c(data.long$pcrc1age_trans, rep(0,n5+n2)), 
                               l.Cumgap = c(data.long$cumgap, rep(0,n5+n2)), 
                               z.Time = c(rep(0,n4), log(data.long.z$gap), rep(0,n2)), # related to log(mu): negative binomial component for longitudinal outcomes
                               intercept.z = c(rep(0,n4), rep(1,n5), rep(0,n2)),                               
                               z.Sex = c(rep(0,n4), data.long.z$sex01, rep(0,n2)),
                               z.Age = c(rep(0,n4), data.long.z$age00_trans, rep(0,n2)),
                               z.Pcrc1age = c(rep(0,n4), data.long.z$pcrc1age_trans, rep(0,n2)),
                               z.Cumgap = c(rep(0,n4), data.long.z$cumgap, rep(0,n2)),
                               d.Sex = c(rep(0,n4+n5), data.term$sex01), # related to hazard function for terminal event
                               d.Adenoma = c(rep(0,n4+n5), data.term$adenoma00),
                               d.Age = c(rep(0,n4+n5), data.term$age00_trans),
                               d.Pcrc1age = c(rep(0,n4+n5), data.term$pcrc1age_trans))

random.covariate <- list(l.YL = c(data.long$id, rep(NA,n5+n2)), 
                         l.ZL = c(data.long$famid, rep(NA,n5+n2)), 
                         z.YL = c(rep(NA,n4), data.long.z$id, rep(NA,n2)), 
                         z.ZL = c(rep(NA,n4), data.long.z$famid, rep(NA,n2)), 
                         d.YL = c(rep(NA,n4+n5), (1:n2)), 
                         d.ZL = c(rep(NA,n4+n5), data.term$famid))

data.f <- c(linear.covariate,random.covariate)
data.f$Y <- y.joint

### Bivariate joint model of longitudinal outcomes and the terminal event
formula1 = Y ~ -1 + mu + intercept.l + l.Sex + l.Age + l.Pcrc1age + l.Cumgap + offset(l.Time) +
  intercept.z + z.Sex + z.Age + z.Pcrc1age + z.Cumgap + offset(z.Time) +
  d.Sex + d.Adenoma + d.Age + d.Pcrc1age +
  f(d.YL, model="iid", n=n2, hyper=list(prec = list(param = c(2.5, 0.5)))) +
  f(l.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(d.ZL, model="iid", n=n3, hyper=list(prec = list(param = c(2.5, 0.5)))) +
  f(l.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) 

inla.model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,control.fixed = list(prec.intercept = 0.1))

### Final Results
summary(inla.model1)
#Time used:
Pre = 1.42, Running = 1.65, Post = 0.0388, Total = 3.1 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu          -2.873 0.572     -3.995   -2.873     -1.752 -2.873   0
intercept.l  1.299 0.430      0.456    1.299      2.142  1.299   0
l.Sex       -0.068 0.418     -0.888   -0.068      0.753 -0.068   0
l.Age       -0.124 0.184     -0.484   -0.124      0.236 -0.124   0
l.Pcrc1age  -0.095 0.618     -1.307   -0.095      1.117 -0.095   0
l.Cumgap    -0.127 0.043     -0.212   -0.127     -0.043 -0.127   0
intercept.z -1.497 0.608     -2.689   -1.497     -0.306 -1.497   0
z.Sex        1.157 0.560      0.060    1.157      2.254  1.157   0
z.Age        0.535 0.259      0.026    0.535      1.043  0.535   0
z.Pcrc1age  -0.879 0.671     -2.194   -0.879      0.437 -0.879   0
z.Cumgap    -0.263 0.082     -0.425   -0.263     -0.102 -0.263   0
d.Sex        2.784 0.664      1.483    2.784      4.085  2.784   0
d.Adenoma   -1.129 0.852     -2.800   -1.129      0.541 -1.129   0
d.Age        1.243 0.270      0.714    1.243      1.773  1.243   0
d.Pcrc1age  -0.392 0.763     -1.887   -0.392      1.103 -0.392   0

Random effects:
  Name	  Model
d.YL IID model
d.ZL IID model
l.YL Copy
z.YL Copy
l.ZL Copy
z.ZL Copy

Model hyperparameters:
  mean       sd 0.025quant 0.5quant 0.975quant   mode
size for nbinomial_0 zero-inflated observations[2] 394.637 4186.692      1.524   22.975   2381.626  3.006
alpha parameter for weibullsurv[3]                   0.494    0.053      0.399    0.491      0.608  0.483
Precision for d.YL                                   1.469    1.750      0.136    0.941      6.021  0.364
Precision for d.ZL                                   3.872    2.314      1.133    3.321      9.877  2.451
Beta for l.YL                                       -0.952    0.483     -1.888   -0.957      0.012 -0.977
Beta for z.YL                                        0.976    0.506     -0.028    0.979      1.966  0.989
Beta for l.ZL                                        1.606    0.531      0.566    1.604      2.659  1.595
Beta for z.ZL                                        1.097    0.533      0.064    1.092      2.164  1.067

Deviance Information Criterion (DIC) ...............: 480.51
Deviance Information Criterion (DIC, saturated) ....: NA
Effective number of parameters .....................: 60.49

Watanabe-Akaike information criterion (WAIC) ...: 492.71
Effective number of parameters .................: 63.87

Marginal log-Likelihood:  -325.73 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model1$marginals.fixed)
#[1] "mu"          "intercept.l" "l.Sex"       "l.Age"       "l.Pcrc1age"  "l.Cumgap"    "intercept.z" "z.Sex"       "z.Age"       "z.Pcrc1age" 
[11] "z.Cumgap"    "d.Sex"       "d.Adenoma"   "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[1]]))
#Mean            20.7886 
Stdev           12.6904 
Quantile  0.025 5.76378 
Quantile  0.25  12.0063 
Quantile  0.5   17.6649 
Quantile  0.75  25.9841 
Quantile  0.975 54.0118 

names(inla.model1$internal.marginals.hyperpar) 
#[1] "log size for nbinomial zero-inflated observations[2]" "alpha_intern for weibullsurv[3]"                     
[3] "Log precision for d.YL"                               "Log precision for d.ZL"                              
[5] "Beta_intern for l.YL"                                 "Beta_intern for z.YL"                                
[7] "Beta_intern for l.ZL"                                 "Beta_intern for z.ZL"     

# for dispersion
inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model1$internal.marginals.hyperpar[[1]]))
#Mean            323.758 
Stdev           2077.77 
Quantile  0.025 -105.935 
Quantile  0.25  0.366493 
Quantile  0.5   34.7979 
Quantile  0.75  171.149 
Quantile  0.975 3445.84 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[3]]))
#Mean            1.16229 
Stdev           0.593199 
Quantile  0.025 0.409371 
Quantile  0.25  0.745839 
Quantile  0.5   1.02821 
Quantile  0.75  1.42786 
Quantile  0.975 2.68618 

# for family-specific random effects Z
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[4]]))
#Mean            0.569082 
Stdev           0.15775 
Quantile  0.025 0.319467 
Quantile  0.25  0.45546 
Quantile  0.5   0.548579 
Quantile  0.75  0.659853 
Quantile  0.975 0.93511 

inla.model1$dic$dic
#[1]  480.5105
inla.model1$dic$family.dic
#[1]  334.396473 142.583516   3.530523
