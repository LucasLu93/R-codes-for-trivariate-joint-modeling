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
                               z.Time = c(rep(0,n4), log(data.long.z$gap), rep(0,n1+n2)), # related to log(mu): negative binomial component for longitudinal outcomes
                               intercept.z = c(rep(0,n4), rep(1,n5), rep(0,n1+n2)),                               
                               z.Sex = c(rep(0,n4), data.long.z$sex01, rep(0,n1+n2)),
                               z.Age = c(rep(0,n4), data.long.z$age00_trans, rep(0,n1+n2)),
                               z.Pcrc1age = c(rep(0,n4), data.long.z$pcrc1age_trans, rep(0,n1+n2)),
                               z.Cumgap = c(rep(0,n4), data.long.z$cumgap, rep(0,n1+n2)),
                               r.Sex = c(rep(0,n4+n5), data$sex01, rep(0,n2)), # relative to hazard function for recurrent events
                               r.Adenoma = c(rep(0,n4+n5), data$adenoma0, rep(0,n2)),
                               r.Age = c(rep(0,n4+n5), data$age0_trans, rep(0,n2)),
                               r.Pcrc1age = c(rep(0,n4+n5), data$pcrc1age_trans, rep(0,n2)),
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

inla.model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,control.fixed = list(prec.intercept = 0.1))

### Final Results
summary(inla.model1)
#Time used:
Pre = 2.26, Running = 8.69, Post = 0.0909, Total = 11 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu1          1.460 0.152      1.163    1.460      1.757  1.460   0
mu2         -3.344 0.585     -4.491   -3.344     -2.198 -3.344   0
intercept.l  1.621 0.413      0.812    1.621      2.430  1.621   0
l.Sex       -0.319 0.408     -1.118   -0.319      0.480 -0.319   0
l.Age       -0.131 0.177     -0.477   -0.131      0.215 -0.131   0
l.Pcrc1age   0.022 0.590     -1.134    0.022      1.178  0.022   0
l.Cumgap    -0.154 0.042     -0.236   -0.154     -0.072 -0.154   0
intercept.z -1.932 0.604     -3.115   -1.932     -0.749 -1.932   0
z.Sex        1.588 0.537      0.536    1.588      2.640  1.588   0
z.Age        0.628 0.246      0.147    0.628      1.110  0.628   0
z.Pcrc1age  -1.075 0.686     -2.419   -1.075      0.269 -1.075   0
z.Cumgap    -0.200 0.077     -0.351   -0.200     -0.049 -0.200   0
r.Sex        0.274 0.149     -0.018    0.274      0.566  0.274   0
r.Adenoma    0.075 0.153     -0.226    0.075      0.375  0.075   0
r.Age        0.008 0.061     -0.112    0.008      0.128  0.008   0
r.Pcrc1age  -0.296 0.264     -0.815   -0.296      0.222 -0.296   0
d.Sex        3.470 0.737      2.026    3.470      4.914  3.470   0
d.Adenoma   -3.008 0.862     -4.699   -3.008     -1.318 -3.008   0
d.Age        1.376 0.317      0.755    1.376      1.997  1.376   0
d.Pcrc1age  -1.090 0.888     -2.831   -1.090      0.650 -1.090   0

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
  mean     sd 0.025quant 0.5quant 0.975quant   mode
size for nbinomial_0 zero-inflated observations[2]  8.625 10.577      1.032    5.487     35.843  2.533
alpha parameter for weibullsurv[3]                  1.192  0.064      1.071    1.191      1.323  1.188
alpha parameter for weibullsurv[4]                  0.699  0.070      0.580    0.693      0.855  0.674
Theta1 for d.YR                                     0.543  0.406     -0.258    0.543      1.339  0.547
Theta2 for d.YR                                    -0.547  0.334     -1.197   -0.549      0.118 -0.560
Theta3 for d.YR                                    -1.267  0.487     -2.244   -1.261     -0.324 -1.237
Theta1 for r.ZR                                     0.930  0.203      0.530    0.931      1.328  0.932
Theta2 for r.ZR                                     0.384  0.306     -0.213    0.383      0.990  0.377
Theta3 for r.ZR                                    -0.757  0.504     -1.767   -0.750      0.217 -0.724
Beta for r.YR                                       0.376  0.160      0.058    0.377      0.689  0.380
Beta for l.YL                                      -0.585  0.258     -1.098   -0.583     -0.083 -0.575
Beta for z.YL                                       0.660  0.285      0.102    0.659      1.223  0.654
Beta for l.ZL                                       1.222  0.450      0.351    1.217      2.122  1.196
Beta for z.ZL                                       1.167  0.430      0.338    1.162      2.030  1.138
Beta for d.ZR                                      -0.532  0.610     -1.758   -0.524      0.644 -0.490

Deviance Information Criterion (DIC) ...............: 146.75
Deviance Information Criterion (DIC, saturated) ....: NA
Effective number of parameters .....................: 143.58

Watanabe-Akaike information criterion (WAIC) ...: 222.50
Effective number of parameters .................: 169.86

Marginal log-Likelihood:  -220.04 
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
#Mean            0.234842 
Stdev           0.0355696 
Quantile  0.025 0.172654 
Quantile  0.25  0.209585 
Quantile  0.5   0.232132 
Quantile  0.75  0.257103 
Quantile  0.975 0.312094

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[2]]))
#Mean            33.5388 
Stdev           20.9957 
Quantile  0.025 9.00176 
Quantile  0.25  19.0635 
Quantile  0.5   28.2916 
Quantile  0.75  41.9755 
Quantile  0.975 88.6871 

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
#Mean            8.53913 
Stdev           9.8268 
Quantile  0.025 1.01191 
Quantile  0.25  2.99007 
Quantile  0.5   5.41267 
Quantile  0.75  10.1667 
Quantile  0.975 35.1169 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model1, "d.YR", return.cov=FALSE) 
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Y.mean <- matrix(rowMeans(sdcor.Y),2,2)
round(sdcor.Y.mean,3)
#[,1]   [,2]
[1,] 1.458 0.875
[2,] 0.875 1.804

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.573 0.122
[2,] 0.122 0.618

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,] 0.546 0.605
[2,] 0.605 0.893

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 2.803 0.970
[2,] 0.970 3.281

# for family-specific random effects Z
mcsamples <- inla.iidkd.sample(10^4, inla.model1, "r.ZR", return.cov=FALSE) 
sdcor.Z <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Z.mean <- matrix(rowMeans(sdcor.Z),2,2)
round(sdcor.Z.mean,3)
#[,1]  [,2]
[1,] 0.461 0.390
[2,] 0.390 0.713

# sd
sdcor.Z.sd <- matrix(apply(sdcor.Z, 1, sd),2,2)
round(sdcor.Z.sd,3)
#[,1]  [,2]
[1,] 0.102 0.258
[2,] 0.258 0.225

# 2.5% lower CI
sdcor.Z.lowerci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Z.lowerci,3)
#[,1]   [,2]
[1,]  0.301 -0.198
[2,] -0.198  0.367

# 97.5% upper CI
sdcor.Z.upperci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Z.upperci,3)
#[,1]  [,2]
[1,] 0.695 0.791
[2,] 0.791 1.244

inla.model1$dic$dic
#[1]  146.7533
inla.model1$dic$family.dic
#[1]  341.12830  144.54031 -287.41019  -51.50509

################################################################################
# zero-inflated Poisson distribution for longitudinal outcomes
inla.model2 <- inla(formula1, family = c("binomial","zeroinflatedpoisson0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,control.fixed = list(prec.intercept = 0.1))

### Final Results
summary(inla.model2)
#Time used:
Pre = 2.25, Running = 7.97, Post = 0.122, Total = 10.3 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu1          1.459 0.152      1.162    1.459      1.756  1.459   0
mu2         -3.298 0.577     -4.429   -3.298     -2.166 -3.298   0
intercept.l  1.614 0.414      0.802    1.614      2.425  1.614   0
l.Sex       -0.320 0.408     -1.119   -0.320      0.479 -0.320   0
l.Age       -0.131 0.177     -0.478   -0.131      0.215 -0.131   0
l.Pcrc1age   0.019 0.593     -1.143    0.019      1.181  0.019   0
l.Cumgap    -0.154 0.042     -0.236   -0.154     -0.072 -0.154   0
intercept.z -1.778 0.559     -2.875   -1.778     -0.682 -1.778   0
z.Sex        1.580 0.499      0.602    1.580      2.558  1.580   0
z.Age        0.627 0.232      0.173    0.627      1.081  0.627   0
z.Pcrc1age  -1.018 0.631     -2.255   -1.018      0.219 -1.018   0
z.Cumgap    -0.222 0.074     -0.367   -0.222     -0.077 -0.222   0
r.Sex        0.277 0.149     -0.016    0.277      0.570  0.277   0
r.Adenoma    0.073 0.153     -0.228    0.073      0.374  0.073   0
r.Age        0.008 0.061     -0.113    0.008      0.128  0.008   0
r.Pcrc1age  -0.299 0.265     -0.817   -0.299      0.220 -0.299   0
d.Sex        3.436 0.729      2.007    3.436      4.866  3.436   0
d.Adenoma   -2.965 0.835     -4.602   -2.965     -1.328 -2.965   0
d.Age        1.356 0.314      0.741    1.356      1.970  1.356   0
d.Pcrc1age  -1.080 0.877     -2.800   -1.080      0.639 -1.080   0

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
  mean    sd 0.025quant 0.5quant 0.975quant   mode
alpha parameter for weibullsurv[3]  1.194 0.065      1.071    1.192      1.328  1.188
alpha parameter for weibullsurv[4]  0.703 0.045      0.622    0.701      0.797  0.695
Theta1 for d.YR                     0.419 0.567     -0.727    0.430      1.505  0.474
Theta2 for d.YR                    -0.514 0.354     -1.190   -0.520      0.201 -0.548
Theta3 for d.YR                    -1.201 0.449     -2.097   -1.196     -0.329 -1.179
Theta1 for r.ZR                     0.914 0.192      0.532    0.915      1.286  0.921
Theta2 for r.ZR                     0.411 0.293     -0.162    0.410      0.990  0.405
Theta3 for r.ZR                    -0.686 0.414     -1.498   -0.687      0.132 -0.690
Beta for r.YR                       0.354 0.168      0.014    0.357      0.676  0.370
Beta for l.YL                      -0.599 0.261     -1.123   -0.596     -0.096 -0.581
Beta for z.YL                       0.686 0.307      0.095    0.682      1.302  0.664
Beta for l.ZL                       1.258 0.421      0.437    1.255      2.097  1.242
Beta for z.ZL                       1.076 0.435      0.219    1.076      1.932  1.077
Beta for d.ZR                      -0.451 0.578     -1.605   -0.446      0.672 -0.424

Deviance Information Criterion (DIC) ...............: 147.63
Deviance Information Criterion (DIC, saturated) ....: 1226.57
Effective number of parameters .....................: 146.89

Watanabe-Akaike information criterion (WAIC) ...: 233.07
Effective number of parameters .................: 178.07

Marginal log-Likelihood:  -220.65 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model2$marginals.fixed)
# [1] "mu1"         "mu2"         "intercept.l" "l.Sex"       "l.Age"       "l.Pcrc1age"  "l.Cumgap"    "intercept.z" "z.Sex"       "z.Age"       "z.Pcrc1age" 
[12] "z.Cumgap"    "r.Sex"       "r.Adenoma"   "r.Age"       "r.Pcrc1age"  "d.Sex"       "d.Adenoma"   "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model2$marginals.fixed[[1]]))
#Mean            0.235094 
Stdev           0.0356024 
Quantile  0.025 0.172847 
Quantile  0.25  0.209814 
Quantile  0.5   0.232381 
Quantile  0.75  0.257375 
Quantile  0.975 0.312416 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model2$marginals.fixed[[2]]))
#Mean            31.8778 
Stdev           19.6662 
Quantile  0.025 8.71964 
Quantile  0.25  18.289 
Quantile  0.5   27.0059 
Quantile  0.75  39.8672 
Quantile  0.975 83.4355 

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
[1,] 1.613 0.852
[2,] 0.852 1.754

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.796 0.13
[2,] 0.130 0.63

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,] 0.527 0.551
[2,] 0.551 0.817

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 3.597 0.964
[2,] 0.964 3.254

# for family-specific random effects Z
mcsamples <- inla.iidkd.sample(10^4, inla.model2, "r.ZR", return.cov=FALSE) 
sdcor.Z <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Z.mean <- matrix(rowMeans(sdcor.Z),2,2)
round(sdcor.Z.mean,3)
#[,1]  [,2]
[1,] 0.459 0.371
[2,] 0.371 0.689

# sd
sdcor.Z.sd <- matrix(apply(sdcor.Z, 1, sd),2,2)
round(sdcor.Z.sd,3)
#[,1]  [,2]
[1,] 0.091 0.211
[2,] 0.211 0.210

# 2.5% lower CI
sdcor.Z.lowerci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Z.lowerci,3)
#[,1]   [,2]
[1,]  0.311 -0.110
[2,] -0.110  0.371

# 97.5% upper CI
sdcor.Z.upperci <- matrix(apply(sdcor.Z,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Z.upperci,3)
#[,1]  [,2]
[1,] 0.665 0.721
[2,] 0.721 1.178

inla.model2$dic$dic
#[1]  147.6263
inla.model2$dic$family.dic
#[1]  341.22311  146.39226 -287.36386  -52.62517

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
  f(d.YR, model="iidkd", order=2, n=2*n2, hyper = list(theta1 = list(param = c(5, 1, 1, 0)))) +
  f(d.YL, copy="d.YR") +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(l.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model3 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,control.fixed = list(prec.intercept = 0.1))

### Final Results
summary(inla.model3)
# Time used:
Pre = 1.58, Running = 7.78, Post = 0.0673, Total = 9.43 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu1          1.301 0.097      1.111    1.301      1.491  1.301   0
mu2         -3.097 0.522     -4.120   -3.097     -2.074 -3.097   0
intercept.l  1.348 0.365      0.633    1.348      2.064  1.348   0
l.Sex       -0.200 0.432     -1.046   -0.200      0.647 -0.200   0
l.Age       -0.261 0.188     -0.629   -0.261      0.108 -0.261   0
l.Pcrc1age   0.040 0.442     -0.826    0.040      0.907  0.040   0
l.Cumgap    -0.146 0.043     -0.231   -0.146     -0.061 -0.146   0
intercept.z -1.375 0.526     -2.406   -1.375     -0.343 -1.375   0
z.Sex        1.200 0.513      0.195    1.200      2.205  1.200   0
z.Age        0.362 0.238     -0.105    0.362      0.828  0.362   0
z.Pcrc1age  -0.658 0.553     -1.742   -0.658      0.426 -0.658   0
z.Cumgap    -0.315 0.079     -0.470   -0.315     -0.160 -0.315   0
r.Sex        0.375 0.162      0.057    0.375      0.693  0.375   0
r.Adenoma   -0.004 0.161     -0.320   -0.004      0.312 -0.004   0
r.Age        0.068 0.063     -0.054    0.068      0.191  0.068   0
r.Pcrc1age  -0.484 0.176     -0.829   -0.484     -0.139 -0.484   0
d.Sex        3.262 0.694      1.902    3.262      4.622  3.262   0
d.Adenoma   -2.775 0.838     -4.416   -2.775     -1.133 -2.775   0
d.Age        1.284 0.296      0.705    1.284      1.864  1.284   0
d.Pcrc1age  -1.124 0.796     -2.684   -1.124      0.435 -1.124   0

Random effects:
  Name	  Model
d.YR IIDKD model
d.YL Copy
r.YR Copy
l.YL Copy
z.YL Copy

Model hyperparameters:
  mean      sd 0.025quant 0.5quant 0.975quant   mode
size for nbinomial_0 zero-inflated observations[2] 35.779 158.908      0.541    8.132    239.738  1.199
alpha parameter for weibullsurv[3]                  1.175   0.053      1.075    1.173      1.283  1.170
alpha parameter for weibullsurv[4]                  0.709   0.037      0.642    0.707      0.787  0.702
Theta1 for d.YR                                    -0.789   0.394     -1.532   -0.799      0.021 -0.848
Theta2 for d.YR                                    -0.115   0.490     -1.107   -0.107      0.823 -0.068
Theta3 for d.YR                                    -0.297   0.295     -0.899   -0.290      0.261 -0.258
Beta for r.YR                                       0.290   0.092      0.117    0.288      0.477  0.277
Beta for l.YL                                      -1.170   0.606     -2.354   -1.174      0.033 -1.188
Beta for z.YL                                       0.995   0.441      0.133    0.993      1.870  0.984

Deviance Information Criterion (DIC) ...............: 191.13
Deviance Information Criterion (DIC, saturated) ....: NA
Effective number of parameters .....................: 167.50

Watanabe-Akaike information criterion (WAIC) ...: 277.57
Effective number of parameters .................: 198.27

Marginal log-Likelihood:  -240.69 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model3$marginals.fixed)
# [1] "mu1"         "mu2"         "intercept.l" "l.Sex"       "l.Age"       "l.Pcrc1age"  "l.Cumgap"    "intercept.z" "z.Sex"       "z.Age"       "z.Pcrc1age" 
[12] "z.Cumgap"    "r.Sex"       "r.Adenoma"   "r.Age"       "r.Pcrc1age"  "d.Sex"       "d.Adenoma"   "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model3$marginals.fixed[[1]]))
#Mean            0.273593 
Stdev           0.0264476 
Quantile  0.025 0.225284 
Quantile  0.25  0.255042 
Quantile  0.5   0.272279 
Quantile  0.75  0.290681 
Quantile  0.975 0.329078 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model3$marginals.fixed[[2]]))
#Mean            25.3141 
Stdev           13.9331 
Quantile  0.025 7.9617 
Quantile  0.25  15.5428 
Quantile  0.5   22.1042 
Quantile  0.75  31.4302 
Quantile  0.975 61.2693 

names(inla.model3$internal.marginals.hyperpar) 
#[1] "log size for nbinomial zero-inflated observations[2]" "alpha_intern for weibullsurv[3]"                     
[3] "alpha_intern for weibullsurv[4]"                      "Theta1 for d.YR"                                     
[5] "Theta2 for d.YR"                                      "Theta3 for d.YR"                                     
[7] "Beta_intern for r.YR"                                 "Beta_intern for l.YL"                                
[9] "Beta_intern for z.YL"    

# for dispersion
inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model3$internal.marginals.hyperpar[[1]]))
#Mean            33.3984 
Stdev           104.841 
Quantile  0.025 0.163178 
Quantile  0.25  2.50448 
Quantile  0.5   7.62458 
Quantile  0.75  24.2314 
Quantile  0.975 235.127 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
mcsamples <- inla.iidkd.sample(10^4, inla.model3, "d.YR", return.cov=FALSE) 
sdcor.Y <- matrix(unlist(mcsamples), nrow = 2^2)
# mean
sdcor.Y.mean <- matrix(rowMeans(sdcor.Y),2,2)
round(sdcor.Y.mean,3)
#[,1]   [,2]
[1,] 2.581 0.306
[2,] 0.306 1.287

# sd
sdcor.Y.sd <- matrix(apply(sdcor.Y, 1, sd),2,2)
round(sdcor.Y.sd,3)
#[,1]  [,2]
[1,] 0.902 0.303
[2,] 0.303 0.766

# 2.5% lower CI
sdcor.Y.lowerci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.025)),2,2)
round(sdcor.Y.lowerci,3)
#[,1]   [,2]
[1,]  1.278 -0.315
[2,] -0.315  0.411

# 97.5% upper CI
sdcor.Y.upperci <- matrix(apply(sdcor.Y,1,function(x) quantile(x,0.975)),2,2)
round(sdcor.Y.upperci,3)
#[,1]  [,2]
[1,] 4.829 0.863
[2,] 0.863 3.257

inla.model3$dic$dic
# [1]  191.1321
inla.model3$dic$family.dic
# [1]  354.26699  153.19835 -265.70198  -50.63127
