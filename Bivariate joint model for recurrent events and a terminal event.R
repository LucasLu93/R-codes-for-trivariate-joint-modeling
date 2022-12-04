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

y.recu <- inla.surv(time = c(time.R/max(time.R), rep(NA, n2)), event = c(delta.R.fij, rep(NA, n2)))
y.term <- inla.surv(time = c(rep(NA, n1), time.D/max(time.D)), event = c(rep(NA, n1), delta.D.fi))
y.joint <- list(y.recu, y.term)

linear.covariate <- data.frame(mu = as.factor(c(rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                               r.Sex = c(data$sex01, rep(0,n2)), # relative to hazard function for recurrent events
                               r.Adenoma = c(data$adenoma0, rep(0,n2)),
                               r.Age = c(data$age0_trans, rep(0,n2)),
                               r.Pcrc1age = c(data$pcrc1age_trans, rep(0,n2)),
                               intercept.d = c(rep(0,n1), rep(1,n2)),
                               d.Sex = c(rep(0,n1), data.term$sex01), # related to hazard function for terminal event
                               d.Adenoma = c(rep(0,n1), data.term$adenoma00),
                               d.Age = c(rep(0,n1), data.term$age00_trans),
                               d.Pcrc1age = c(rep(0,n1), data.term$pcrc1age_trans))

random.covariate <- list(r.YR = c(data$id, rep(NA,n2)), 
                         r.ZR = c(data$famid, rep(NA,n2)), 
                         d.YR = c(rep(NA,n1), 1:n2), 
                         d.ZR = c(rep(NA,n1), data.term$famid))

data.f <- c(linear.covariate,random.covariate)
data.f$Y <- y.joint

### Bivariate joint model of recurrent events and the terminal event
formula1 = Y ~ -1 + mu + 
  r.Sex + r.Adenoma + r.Age + r.Pcrc1age +
  d.Sex + d.Adenoma + d.Age + d.Pcrc1age +
  f(d.YR, model="iid", n=n2, hyper=list(prec = list(param = c(5, 0.5)))) +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(r.ZR, model="iid", n=n3, hyper=list(prec = list(param = c(5, 0.5)))) +
  f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model1 <- inla(formula1, family = c("weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE)

### Final Results
summary(inla.model1)

#Call:
c("inla.core(formula = formula, family = family, contrasts = contrasts, ", " data = data, quantiles = quantiles, E = E, 
   offset = offset, ", " scale = scale, weights = weights, Ntrials = Ntrials, strata = strata, ", " lp.scale = lp.scale, 
   link.covariates = link.covariates, verbose = verbose, ", " lincomb = lincomb, selection = selection, control.compute = 
   control.compute, ", " control.predictor = control.predictor, control.family = control.family, ", " control.inla = 
   control.inla, control.fixed = control.fixed, ", " control.mode = control.mode, control.expert = control.expert, ", " 
   control.hazard = control.hazard, control.lincomb = control.lincomb, ", " control.update = control.update, control.lp.scale 
   = control.lp.scale, ", " control.pardiso = control.pardiso, only.hyperparam = only.hyperparam, ", " inla.call = inla.call, 
   inla.arg = inla.arg, num.threads = num.threads, ", " blas.num.threads = blas.num.threads, keep = keep, working.directory = 
   working.directory, ", " silent = silent, inla.mode = inla.mode, safe = FALSE, debug = debug, ", " .parent.frame = 
   .parent.frame)") 
Time used:
  Pre = 1.54, Running = 0.458, Post = 0.0192, Total = 2.02 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu1         1.391 0.130      1.135    1.391      1.647  1.391   0
mu2        -2.853 0.538     -3.906   -2.853     -1.799 -2.853   0
r.Sex       0.160 0.139     -0.113    0.160      0.434  0.160   0
r.Adenoma   0.256 0.157     -0.051    0.256      0.563  0.256   0
r.Age       0.008 0.058     -0.106    0.008      0.123  0.008   0
r.Pcrc1age -0.273 0.227     -0.718   -0.273      0.172 -0.273   0
d.Sex       2.752 0.640      1.498    2.752      4.006  2.752   0
d.Adenoma  -0.663 0.832     -2.294   -0.663      0.967 -0.663   0
d.Age       1.124 0.246      0.641    1.124      1.608  1.124   0
d.Pcrc1age -0.270 0.713     -1.668   -0.270      1.127 -0.270   0

Random effects:
  Name	  Model
d.YR IID model
r.ZR IID model
r.YR Copy
d.ZR Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant  mode
alpha parameter for weibullsurv    1.144 0.070      1.012    1.142      1.288 1.137
alpha parameter for weibullsurv[2] 0.490 0.049      0.398    0.488      0.592 0.485
Precision for d.YR                 8.034 3.956      2.665    7.274     17.870 5.847
Precision for r.ZR                 8.341 3.103      3.714    7.865     15.761 6.969
Beta for r.YR                      0.935 0.338      0.236    0.948      1.566 0.998
Beta for d.ZR                      0.088 0.639     -1.148    0.079      1.370 0.045

Deviance Information Criterion (DIC) ...............: 363.43
Deviance Information Criterion (DIC, saturated) ....: 806.88
Effective number of parameters .....................: 363.43

Watanabe-Akaike information criterion (WAIC) ...: -249.79
Effective number of parameters .................: 56.59

Marginal log-Likelihood:  79.76 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model1$marginals.fixed)
#[1] "mu1"        "mu2"        "r.Sex"      "r.Adenoma"  "r.Age"      "r.Pcrc1age" "d.Sex"      "d.Adenoma"  "d.Age"      "d.Pcrc1age"
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[1]]))
#Mean            0.250918 
Stdev           0.0326407 
Quantile  0.025 0.192858 
Quantile  0.25  0.22784 
Quantile  0.5   0.248765 
Quantile  0.75  0.271612 
Quantile  0.975 0.320878 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[2]]))
#Mean            19.9857 
Stdev           11.3717 
Quantile  0.025 6.04441 
Quantile  0.25  12.0416 
Quantile  0.5   17.308 
Quantile  0.75  24.8728 
Quantile  0.975 49.4706 

names(inla.model1$internal.marginals.hyperpar) 
#[1] "alpha_intern for weibullsurv"    "alpha_intern for weibullsurv[2]" "Log precision for d.YR"          "Log precision for r.ZR"         
[5] "Beta_intern for r.YR"            "Beta_intern for d.ZR"      

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[3]]))
#Mean            0.384457 
Stdev           0.0952532 
Quantile  0.025 0.237462 
Quantile  0.25  0.316087 
Quantile  0.5   0.370211 
Quantile  0.75  0.437884 
Quantile  0.975 0.609266 

# for family-specific random effects Z
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[4]]))
#Mean            0.363963 
Stdev           0.067528 
Quantile  0.025 0.252602 
Quantile  0.25  0.315638 
Quantile  0.5   0.356285 
Quantile  0.75  0.404251 
Quantile  0.975 0.516953 

