rm(list=ls())
library(INLA)
inla.setOption(inla.mode="experimental")
inla.setOption(num.threads=1) 

# import the data
data<-read.table("~/Desktop/NFdat_newer.txt", stringsAsFactors=FALSE)

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

l.long <- c(z, rep(NA, n5), rep(NA, n1), rep(NA, n2))
z.long <- c(rep(NA, n4), data.long.z$num, rep(NA, n1), rep(NA, n2))
y.recu <- inla.surv(time = c(rep(NA, n4+n5), time.R/max(time.R), rep(NA, n2)), event = c(rep(NA, n4+n5), delta.R.fij, rep(NA, n2)))
y.term <- inla.surv(time = c(rep(NA, n4+n5+n1), time.D/max(time.D)), event = c(rep(NA, n4+n5+n1), delta.D.fi))
y.joint <- list(l.long, z.long, y.recu, y.term)

linear.covariate <- data.frame(mu = as.factor(c(rep(NA,n4+n5),rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                               l.Time = c(log(data.long$gap), rep(0,n5+n1+n2)), # related to logit(pi): logistic component for longitudinal outcomes 
                               intercept.l = c(rep(1,n4), rep(0,n5+n1+n2)),
                               l.Sex = c(data.long$sex01, rep(0,n5+n1+n2)), 
                               l.Age = c(data.long$age00_trans, rep(0,n5+n1+n2)), 
                               l.Pcrc1age = c(data.long$pcrc1age_trans, rep(0,n5+n1+n2)), 
                               l.Cumgap = c(data.long$cumgap, rep(0,n5+n1+n2)), 
                               l.Polyps = c(data.long$polyps00, rep(0,n5+n1+n2)), 
                               l.Adenoma = c(data.long$adenoma00, rep(0,n5+n1+n2)), 
                               z.Time = c(rep(0,n4), log(data.long.z$gap), rep(0,n1+n2)), # related to log(mu): negative binomial component for longitudinal outcomes
                               intercept.z = c(rep(0,n4), rep(1,n5), rep(0,n1+n2)),                               
                               z.Sex = c(rep(0,n4), data.long.z$sex01, rep(0,n1+n2)),
                               z.Age = c(rep(0,n4), data.long.z$age00_trans, rep(0,n1+n2)),
                               z.Pcrc1age = c(rep(0,n4), data.long.z$pcrc1age_trans, rep(0,n1+n2)),
                               z.Cumgap = c(rep(0,n4), data.long.z$cumgap, rep(0,n1+n2)),
                               z.Polyps = c(rep(0,n4), data.long.z$polyps00, rep(0,n1+n2)),
                               z.Adenoma = c(rep(0,n4), data.long.z$adenoma00, rep(0,n1+n2)),
                               r.Sex = c(rep(0,n4+n5), data$sex01, rep(0,n2)), # relative to hazard function for recurrent events
                               r.Adenoma = c(rep(0,n4+n5), data$adenoma0, rep(0,n2)),
                               r.Age = c(rep(0,n4+n5), data$age0_trans, rep(0,n2)),
                               r.Pcrc1age = c(rep(0,n4+n5), data$pcrc1age_trans, rep(0,n2)),
                               intercept.d = c(rep(0,n4+n5+n1), rep(1,n2)),
                               d.Sex = c(rep(0,n4+n5+n1), data.term$sex01), # related to hazard function for terminal event
                               d.Adenoma = c(rep(0,n4+n5+n1), data.term$adenoma00),
                               d.Age = c(rep(0,n4+n5+n1), data.term$age00_trans),
                               d.Pcrc1age = c(rep(0,n4+n5+n1), data.term$pcrc1age_trans))

random.covariate <- list(l.YL = c(data.long$id+n2, rep(NA,n5+n1+n2)), 
                         l.ZL = c(data.long$famid+n3, rep(NA,n5+n1+n2)), 
                         z.YL = c(rep(NA,n4), data.long.z$id+n2, rep(NA,n1+n2)), 
                         z.ZL = c(rep(NA,n4), data.long.z$famid+n3, rep(NA,n1+n2)), 
                         r.YR = c(rep(NA,n4+n5), data$id, rep(NA,n2)), 
                         r.ZR = c(rep(NA,n4+n5), data$famid, rep(NA,n2)), 
                         d.YR = c(rep(NA,n4+n5+n1), 1:n2), 
                         d.ZR = c(rep(NA,n4+n5+n1), data.term$famid),
                         d.YL = c(rep(NA,n4+n5+n1), (1:n2)+n2), 
                         d.ZL = c(rep(NA,n4+n5+n1), data.term$famid+n3))

data.f <- c(linear.covariate,random.covariate)
data.f$Y <- y.joint

### Trivariate joint model
formula1 = Y ~ -1 + mu + intercept.l + l.Sex + l.Age + l.Pcrc1age + l.Cumgap + offset(l.Time) +
  intercept.z + z.Sex + z.Age + z.Pcrc1age + z.Cumgap + offset(z.Time) +
  r.Sex + r.Adenoma + r.Age + r.Pcrc1age +
  d.Sex + d.Adenoma + d.Age + d.Pcrc1age +
  f(d.YR, model="iidkd", order=2, n=2*n2, hyper = list(theta1 = list(param = c(10, 1, 1, 0)))) +
  f(d.YL, copy="d.YR") +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(l.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(r.ZR, model="iidkd", order=2, n=2*n3, hyper = list(theta1 = list(param = c(10, 1, 1, 0)))) +
  f(d.ZL, copy="r.ZR") +
  f(l.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.ZL, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE)

### Final Results
summary(inla.model1)
#Call:
c("inla.core(formula = formula, family = family, contrasts = contrasts, ", " data = data, quantiles = 
   quantiles, E = E, offset = offset, ", " scale = scale, weights = weights, Ntrials = Ntrials, strata = 
   strata, ", " lp.scale = lp.scale, link.covariates = link.covariates, verbose = verbose, ", " lincomb = 
   lincomb, selection = selection, control.compute = control.compute, ", " control.predictor = 
   control.predictor, control.family = control.family, ", " control.inla = control.inla, control.fixed = 
   control.fixed, ", " control.mode = control.mode, control.expert = control.expert, ", " control.hazard = 
   control.hazard, control.lincomb = control.lincomb, ", " control.update = control.update, 
   control.lp.scale = control.lp.scale, ", " control.pardiso = control.pardiso, only.hyperparam = 
   only.hyperparam, ", " inla.call = inla.call, inla.arg = inla.arg, num.threads = num.threads, ", " 
   blas.num.threads = blas.num.threads, keep = keep, working.directory = working.directory, ", " silent = 
   silent, inla.mode = inla.mode, safe = FALSE, debug = debug, ", " .parent.frame = .parent.frame)") 
Time used:
  Pre = 4.93, Running = 7.27, Post = 0.0597, Total = 12.3 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant mode kld
mu1          1.420 0.133      1.159    1.420      1.681   NA   0
mu2         -2.713 0.524     -3.740   -2.713     -1.685   NA   0
intercept.l  1.457 0.400      0.673    1.457      2.240   NA   0
l.Sex       -0.134 0.394     -0.907   -0.134      0.639   NA   0
l.Age       -0.138 0.174     -0.480   -0.138      0.203   NA   0
l.Pcrc1age   0.002 0.559     -1.093    0.002      1.097   NA   0
l.Cumgap    -0.136 0.042     -0.218   -0.136     -0.054   NA   0
intercept.z -1.541 0.567     -2.653   -1.541     -0.429   NA   0
z.Sex        1.100 0.525      0.070    1.100      2.129   NA   0
z.Age        0.452 0.244     -0.026    0.452      0.930   NA   0
z.Pcrc1age  -0.755 0.605     -1.941   -0.755      0.431   NA   0
z.Cumgap    -0.243 0.079     -0.398   -0.243     -0.089   NA   0
r.Sex        0.187 0.142     -0.092    0.187      0.466   NA   0
r.Adenoma    0.150 0.155     -0.154    0.150      0.455   NA   0
r.Age        0.010 0.059     -0.106    0.010      0.126   NA   0
r.Pcrc1age  -0.273 0.232     -0.728   -0.273      0.182   NA   0
d.Sex        2.799 0.629      1.567    2.799      4.031   NA   0
d.Adenoma   -1.472 0.808     -3.055   -1.472      0.111   NA   0
d.Age        1.192 0.255      0.691    1.192      1.692   NA   0
d.Pcrc1age  -0.554 0.716     -1.958   -0.554      0.849   NA   0

Random effects:
  Name	  Model
d.YR IIDKD model
r.ZR IIDKD model
d.YL Copy
r.YR Copy
l.YL Copy
z.YL Copy
d.ZL Copy
l.ZL Copy
z.ZL Copy
d.ZR Copy

Model hyperparameters:
  mean     sd 0.025quant 0.5quant 0.975quant mode
size for nbinomial zero-inflated observations[2] 14.017 10.041      3.161   11.396     40.992   NA
alpha parameter for weibullsurv[3]                1.167  0.057      1.059    1.165      1.282   NA
alpha parameter for weibullsurv[4]                0.535  0.051      0.445    0.531      0.647   NA
Theta1 for d.YR                                   0.973  0.245      0.486    0.974      1.454   NA
Theta2 for d.YR                                   0.268  0.299     -0.338    0.271      0.853   NA
Theta3 for d.YR                                  -1.432  0.563     -2.566   -1.424     -0.342   NA
Theta1 for r.ZR                                   1.108  0.157      0.789    1.110      1.415   NA
Theta2 for r.ZR                                   0.787  0.236      0.318    0.789      1.249   NA
Theta3 for r.ZR                                  -1.048  0.523     -2.074   -1.052     -0.003   NA
Beta for r.YR                                     0.750  0.242      0.261    0.753      1.223   NA
Beta for l.YL                                    -1.051  0.403     -1.834   -1.055     -0.243   NA
Beta for z.YL                                     1.283  0.412      0.447    1.287      2.087   NA
Beta for l.ZL                                     1.699  0.450      0.813    1.698      2.588   NA
Beta for z.ZL                                     1.086  0.411      0.287    1.081      1.913   NA
Beta for d.ZR                                    -0.567  0.557     -1.651   -0.572      0.544   NA

Deviance Information Criterion (DIC) ...............: 197.56
Deviance Information Criterion (DIC, saturated) ....: -3888.58
Effective number of parameters .....................: 102.44

Watanabe-Akaike information criterion (WAIC) ...: 229.53
Effective number of parameters .................: 116.66

Marginal log-Likelihood:  -224.43 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model1$marginals.fixed)
#[1] "mu1"         "mu2"         "intercept.l" "l.Sex"       "l.Age"       "l.Pcrc1age"  "l.Cumgap"    "intercept.z"
[9] "z.Sex"       "z.Age"       "z.Pcrc1age"  "z.Cumgap"    "r.Sex"       "r.Adenoma"   "r.Age"       "r.Pcrc1age" 
[17] "d.Sex"       "d.Adenoma"   "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[1]]))
#Mean            0.243877 
Stdev           0.0323746 
Quantile  0.025 0.186413 
Quantile  0.25  0.220975 
Quantile  0.5   0.2417 
Quantile  0.75  0.26437 
Quantile  0.975 0.313383 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[2]]))
#Mean            17.2554 
Stdev           9.54614 
Quantile  0.025 5.3949 
Quantile  0.25  10.5645 
Quantile  0.5   15.0486 
Quantile  0.75  21.4324 
Quantile  0.975 41.908 

names(inla.model1$internal.marginals.hyperpar) 
#[1] "log size for nbinomial zero-inflated observations[2]" "alpha_intern for weibullsurv[3]"                     
[3] "alpha_intern for weibullsurv[4]"                      "Theta1 for d.YR"                                     
[5] "Theta2 for d.YR"                                      "Theta3 for d.YR"                                     
[7] "Theta1 for r.ZR"                                      "Theta2 for r.ZR"                                     
[9] "Theta3 for r.ZR"                                      "Beta_intern for r.YR"                                
[11] "Beta_intern for l.YL"                                 "Beta_intern for z.YL"                                
[13] "Beta_intern for l.ZL"                                 "Beta_intern for z.ZL"                                
[15] "Beta_intern for d.ZR"   

# for dispersion
inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model1$internal.marginals.hyperpar[[1]]))
#Mean            13.9718 
Stdev           9.77217 
Quantile  0.025 3.21483 
Quantile  0.25  7.36064 
Quantile  0.5   11.3739 
Quantile  0.75  17.5608 
Quantile  0.975 39.9843 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model1, "d.YR", return.cov=FALSE) 
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Y.mean <- matrix(rowMeans(sdcor.Y),2,2)
round(sdcor.Y.mean,3)
#[,1]   [,2]
[1,] 0.621 0.676
[2,] 0.676 0.824

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.269 0.206
[2,] 0.206 0.334

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,] 0.280 0.162
[2,] 0.162 0.361

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 1.293 0.934
[2,] 0.934 1.655

# for family-specific random effects Z
mcsamples <- inla.iidkd.sample(10^4, inla.model1, "r.ZR", return.cov=FALSE) 
sdcor.Z <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Z.mean <- matrix(rowMeans(sdcor.Z),2,2)
round(sdcor.Z.mean,3)
#[,1]  [,2]
[1,] 0.387 0.395
[2,] 0.395 0.468

# sd
sdcor.Z.sd <- matrix(apply(sdcor.Z, 1, sd),2,2)
round(sdcor.Z.sd,3)
#[,1]  [,2]
[1,] 0.080 0.226
[2,] 0.226 0.115

# 2.5% lower CI
sdcor.Z.lowerci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Z.lowerci,3)
#[,1]   [,2]
[1,]  0.260 -0.108
[2,] -0.108  0.283

# 97.5% upper CI
sdcor.Z.upperci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Z.upperci,3)
#[,1]  [,2]
[1,] 0.568 0.759
[2,] 0.759 0.727

################################################################################
# zero-inflated Poisson distribution for longitudinal outcomes
inla.model2 <- inla(formula1, family = c("binomial","zeroinflatedpoisson0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE)

### Final Results
summary(inla.model2)
#Call:
c("inla.core(formula = formula, family = family, contrasts = contrasts, ", " data = data, quantiles = quantiles, E = E, offset = offset, ", " scale = 
   scale, weights = weights, Ntrials = Ntrials, strata = strata, ", " lp.scale = lp.scale, link.covariates = link.covariates, verbose = verbose, ", " lincomb 
   = lincomb, selection = selection, control.compute = control.compute, ", " control.predictor = control.predictor, control.family = control.family, ", " 
   control.inla = control.inla, control.fixed = control.fixed, ", " control.mode = control.mode, control.expert = control.expert, ", " control.hazard = 
   control.hazard, control.lincomb = control.lincomb, ", " control.update = control.update, control.lp.scale = control.lp.scale, ", " control.pardiso = 
   control.pardiso, only.hyperparam = only.hyperparam, ", " inla.call = inla.call, inla.arg = inla.arg, num.threads = num.threads, ", " blas.num.threads = 
   blas.num.threads, keep = keep, working.directory = working.directory, ", " silent = silent, inla.mode = inla.mode, safe = FALSE, debug = debug, ", " 
   .parent.frame = .parent.frame)") 
Time used:
  Pre = 5.37, Running = 6.59, Post = 0.0656, Total = 12 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant mode kld
mu1          1.421 0.133      1.160    1.421      1.681   NA   0
mu2         -2.686 0.530     -3.724   -2.686     -1.648   NA   0
intercept.l  1.476 0.400      0.693    1.476      2.260   NA   0
l.Sex       -0.146 0.396     -0.923   -0.146      0.631   NA   0
l.Age       -0.140 0.175     -0.482   -0.140      0.203   NA   0
l.Pcrc1age   0.005 0.559     -1.090    0.005      1.101   NA   0
l.Cumgap    -0.138 0.042     -0.221   -0.138     -0.056   NA   0
intercept.z -1.500 0.548     -2.575   -1.500     -0.425   NA   0
z.Sex        1.102 0.510      0.103    1.102      2.101   NA   0
z.Age        0.455 0.238     -0.012    0.455      0.921   NA   0
z.Pcrc1age  -0.726 0.585     -1.874   -0.726      0.421   NA   0
z.Cumgap    -0.251 0.077     -0.403   -0.251     -0.099   NA   0
r.Sex        0.192 0.143     -0.087    0.192      0.472   NA   0
r.Adenoma    0.142 0.155     -0.162    0.142      0.445   NA   0
r.Age        0.010 0.059     -0.106    0.010      0.126   NA   0
r.Pcrc1age  -0.276 0.232     -0.730   -0.276      0.179   NA   0
d.Sex        2.798 0.633      1.557    2.798      4.039   NA   0
d.Adenoma   -1.546 0.808     -3.130   -1.546      0.037   NA   0
d.Age        1.189 0.258      0.684    1.189      1.694   NA   0
d.Pcrc1age  -0.572 0.719     -1.980   -0.572      0.837   NA   0

Random effects:
  Name	  Model
d.YR IIDKD model
r.ZR IIDKD model
d.YL Copy
r.YR Copy
l.YL Copy
z.YL Copy
d.ZL Copy
l.ZL Copy
z.ZL Copy
d.ZR Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant mode
alpha parameter for weibullsurv[3]  1.167 0.064      1.043    1.166      1.296   NA
alpha parameter for weibullsurv[4]  0.529 0.035      0.460    0.529      0.599   NA
Theta1 for d.YR                     0.949 0.268      0.415    0.951      1.472   NA
Theta2 for d.YR                     0.267 0.347     -0.399    0.260      0.967   NA
Theta3 for d.YR                    -1.401 0.552     -2.499   -1.397     -0.322   NA
Theta1 for r.ZR                     1.110 0.185      0.740    1.112      1.469   NA
Theta2 for r.ZR                     0.795 0.236      0.326    0.796      1.255   NA
Theta3 for r.ZR                    -1.066 0.571     -2.205   -1.062      0.051   NA
Beta for r.YR                       0.744 0.290      0.171    0.744      1.315   NA
Beta for l.YL                      -1.079 0.480     -2.042   -1.073     -0.146   NA
Beta for z.YL                       1.292 0.416      0.491    1.284      2.132   NA
Beta for l.ZL                       1.743 0.501      0.779    1.734      2.754   NA
Beta for z.ZL                       1.042 0.534     -0.011    1.042      2.094   NA
Beta for d.ZR                      -0.609 0.527     -1.664   -0.602      0.413   NA

Deviance Information Criterion (DIC) ...............: 196.78
Deviance Information Criterion (DIC, saturated) ....: -3530.38
Effective number of parameters .....................: 105.67

Watanabe-Akaike information criterion (WAIC) ...: 233.01
Effective number of parameters .................: 122.29

Marginal log-Likelihood:  -225.24 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model2$marginals.fixed)
#[1] "mu1"         "mu2"         "intercept.l" "l.Sex"       "l.Polyps"    "l.Age"       "l.Pcrc1age"  "intercept.z" "z.Sex"       "z.Polyps"    "z.Age"      
[12] "z.Pcrc1age"  "r.Sex"       "r.Polyps"    "r.Age"       "r.Pcrc1age"  "d.Sex"       "d.Polyps"    "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model2$marginals.fixed[[1]]))
#Mean            0.24363 
Stdev           0.0322826 
Quantile  0.025 0.186318 
Quantile  0.25  0.220794 
Quantile  0.5   0.241463 
Quantile  0.75  0.264067 
Quantile  0.975 0.312928 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model2$marginals.fixed[[2]]))
#Mean            16.8544 
Stdev           9.42856 
Quantile  0.025 5.20128 
Quantile  0.25  10.2544 
Quantile  0.5   14.6588 
Quantile  0.75  20.9511 
Quantile  0.975 41.2422 

names(inla.model2$internal.marginals.hyperpar) 
# [1] "alpha_intern for weibullsurv[3]" "alpha_intern for weibullsurv[4]" "Theta1 for d.YR"                
[4] "Theta2 for d.YR"                 "Theta3 for d.YR"                 "Theta1 for r.ZR"                
[7] "Theta2 for r.ZR"                 "Theta3 for r.ZR"                 "Beta_intern for r.YR"           
[10] "Beta_intern for l.YL"            "Beta_intern for z.YL"            "Beta_intern for l.ZL"           
[13] "Beta_intern for z.ZL"            "Beta_intern for d.ZR"   

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model2, "d.YR", return.cov=FALSE) 
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Y.mean <- matrix(rowMeans(sdcor.Y),2,2)
round(sdcor.Y.mean,3)
#[,1]   [,2]
[1,] 0.620 0.679
[2,] 0.679 0.810

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.244 0.188
[2,] 0.188 0.291

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,] 0.301 0.220
[2,] 0.220 0.376

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 1.237 0.930
[2,] 0.930 1.491

# for family-specific random effects Z
mcsamples <- inla.iidkd.sample(10^4, inla.model2, "r.ZR", return.cov=FALSE) 
sdcor.Z <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Z.mean <- matrix(rowMeans(sdcor.Z),2,2)
round(sdcor.Z.mean,3)
#[,1]  [,2]
[1,] 0.383 0.428
[2,] 0.428 0.471

# sd
sdcor.Z.sd <- matrix(apply(sdcor.Z, 1, sd),2,2)
round(sdcor.Z.sd,3)
#[,1]  [,2]
[1,] 0.075 0.205
[2,] 0.205 0.114

# 2.5% lower CI
sdcor.Z.lowerci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Z.lowerci,3)
#[,1]   [,2]
[1,]  0.264 -0.023
[2,] -0.023  0.287

# 97.5% upper CI
sdcor.Z.upperci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Z.upperci,3)
#[,1]  [,2]
[1,] 0.555 0.765
[2,] 0.765 0.729

################################################################################  
# ignore the familial correlations completely 
random.covariate <- list(l.YL = c(data.long$id+n2, rep(NA,n5+n1+n2)), 
                         z.YL = c(rep(NA,n4), data.long.z$id+n2, rep(NA,n1+n2)), 
                         r.YR = c(rep(NA,n4+n5), data$id, rep(NA,n2)), 
                         d.YR = c(rep(NA,n4+n5+n1), 1:n2), 
                         d.YL = c(rep(NA,n4+n5+n1), (1:n2)+n2))

data.f <- c(linear.covariate,random.covariate)
data.f$Y <- y.joint

### Trivariate joint model
formula1 = Y ~ -1 + mu + intercept.l + l.Sex + l.Age + l.Pcrc1age + l.Cumgap + offset(l.Time) +
  intercept.z + z.Sex + z.Age + z.Pcrc1age + z.Cumgap + offset(z.Time) +
  r.Sex + r.Adenoma + r.Age + r.Pcrc1age +
  d.Sex + d.Adenoma + d.Age + d.Pcrc1age +
  f(d.YR, model="iidkd", order=2, n=2*n2, hyper = list(theta1 = list(param = c(10, 1, 1, 0)))) +
  f(d.YL, copy="d.YR") +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(l.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model3 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE)

### Final Results
summary(inla.model3)
# Call:
c("inla.core(formula = formula, family = family, contrasts = contrasts, ", " data = data, quantiles = quantiles, E = E, 
   offset = offset, ", " scale = scale, weights = weights, Ntrials = Ntrials, strata = strata, ", " lp.scale = lp.scale, 
   link.covariates = link.covariates, verbose = verbose, ", " lincomb = lincomb, selection = selection, control.compute = 
   control.compute, ", " control.predictor = control.predictor, control.family = control.family, ", " control.inla = 
   control.inla, control.fixed = control.fixed, ", " control.mode = control.mode, control.expert = control.expert, ", " 
   control.hazard = control.hazard, control.lincomb = control.lincomb, ", " control.update = control.update, 
   control.lp.scale = control.lp.scale, ", " control.pardiso = control.pardiso, only.hyperparam = only.hyperparam, ", " 
   inla.call = inla.call, inla.arg = inla.arg, num.threads = num.threads, ", " blas.num.threads = blas.num.threads, keep = 
   keep, working.directory = working.directory, ", " silent = silent, inla.mode = inla.mode, safe = FALSE, debug = debug, 
   ", " .parent.frame = .parent.frame)") 
Time used:
  Pre = 3.91, Running = 4.95, Post = 0.0619, Total = 8.92 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant mode kld
mu1          1.246 0.092      1.065    1.246      1.426   NA   0
mu2         -2.953 0.542     -4.016   -2.953     -1.890   NA   0
intercept.l  0.627 0.321     -0.002    0.627      1.257   NA   0
l.Sex        0.254 0.378     -0.486    0.254      0.995   NA   0
l.Age       -0.211 0.161     -0.527   -0.211      0.105   NA   0
l.Pcrc1age  -0.332 0.381     -1.080   -0.332      0.415   NA   0
l.Cumgap    -0.084 0.039     -0.161   -0.084     -0.008   NA   0
intercept.z -0.506 0.574     -1.630   -0.506      0.619   NA   0
z.Sex        0.982 0.557     -0.109    0.982      2.074   NA   0
z.Age        0.396 0.257     -0.108    0.396      0.899   NA   0
z.Pcrc1age  -0.686 0.601     -1.863   -0.686      0.492   NA   0
z.Cumgap    -0.336 0.085     -0.503   -0.336     -0.169   NA   0
r.Sex        0.244 0.156     -0.062    0.244      0.550   NA   0
r.Adenoma    0.208 0.161     -0.109    0.208      0.524   NA   0
r.Age        0.058 0.060     -0.060    0.058      0.176   NA   0
r.Pcrc1age  -0.414 0.165     -0.739   -0.414     -0.090   NA   0
d.Sex        2.928 0.660      1.635    2.928      4.220   NA   0
d.Adenoma   -0.994 0.836     -2.633   -0.994      0.645   NA   0
d.Age        1.157 0.261      0.645    1.157      1.669   NA   0
d.Pcrc1age  -0.450 0.713     -1.847   -0.450      0.948   NA   0

Random effects:
  Name	  Model
d.YR IIDKD model
d.YL Copy
r.YR Copy
l.YL Copy
z.YL Copy

Model hyperparameters:
  mean     sd 0.025quant 0.5quant 0.975quant mode
size for nbinomial zero-inflated observations[2] 33.442 40.072      5.564   29.784    175.718   NA
alpha parameter for weibullsurv[3]                1.136  0.028      1.068    1.131      1.188   NA
alpha parameter for weibullsurv[4]                0.518  0.014      0.492    0.520      0.554   NA
Theta1 for d.YR                                   0.872  0.138      0.548    0.882      1.086   NA
Theta2 for d.YR                                   0.385  0.131      0.172    0.379      0.687   NA
Theta3 for d.YR                                  -1.642  0.289     -2.383   -1.705     -1.097   NA
Beta for r.YR                                     0.939  0.158      0.555    0.923      1.224   NA
Beta for l.YL                                     1.411  0.214      1.049    1.402      1.893   NA
Beta for z.YL                                     1.703  0.158      1.340    1.703      1.983   NA

Deviance Information Criterion (DIC) ...............: 249.56
Deviance Information Criterion (DIC, saturated) ....: -3686.13
Effective number of parameters .....................: 112.39

Watanabe-Akaike information criterion (WAIC) ...: 290.97
Effective number of parameters .................: 130.30

Marginal log-Likelihood:  -245.59 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model3$marginals.fixed)
#[1] "mu1"         "mu2"         "intercept.l" "l.Sex"       "l.Polyps"    "l.Age"       "l.Pcrc1age"  "intercept.z" "z.Sex"       "z.Polyps"    "z.Age"      
[12] "z.Pcrc1age"  "r.Sex"       "r.Polyps"    "r.Age"       "r.Pcrc1age"  "d.Sex"       "d.Polyps"    "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model3$marginals.fixed[[1]]))
#Mean            0.288992 
Stdev           0.0264482 
Quantile  0.025 0.240487 
Quantile  0.25  0.270467 
Quantile  0.5   0.287745 
Quantile  0.75  0.306128 
Quantile  0.975 0.344292 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model3$marginals.fixed[[2]]))
#Mean            22.1576 
Stdev           12.7329 
Quantile  0.025 6.622 
Quantile  0.25  13.2735 
Quantile  0.5   19.1401 
Quantile  0.75  27.5942 
Quantile  0.975 55.2179 

names(inla.model3$internal.marginals.hyperpar) 
#[1] "log size for nbinomial zero-inflated observations[2]" "alpha_intern for weibullsurv[3]"                     
[3] "alpha_intern for weibullsurv[4]"                      "Theta1 for d.YR"                                     
[5] "Theta2 for d.YR"                                      "Theta3 for d.YR"                                     
[7] "Beta_intern for r.YR"                                 "Beta_intern for l.YL"                                
[9] "Beta_intern for z.YL"    

# for dispersion
inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model3$internal.marginals.hyperpar[[1]]))
#Mean            33.0552 
Stdev           36.5973 
Quantile  0.025 7.03644 
Quantile  0.25  13.5705 
Quantile  0.5   21.184 
Quantile  0.75  37.5376 
Quantile  0.975 131.338 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model3, "d.YR", return.cov=FALSE) 
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Y.mean <- matrix(rowMeans(sdcor.Y),2,2)
round(sdcor.Y.mean,3)
#[,1]   [,2]
[1,] 0.668 0.716
[2,] 0.716 0.732

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.201 0.147
[2,] 0.147 0.213

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,] 0.381 0.349
[2,] 0.349 0.401

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 1.150 0.916
[2,] 0.916 1.237

