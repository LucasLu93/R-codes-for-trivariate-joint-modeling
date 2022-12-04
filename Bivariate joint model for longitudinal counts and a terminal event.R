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
                               l.Polyps = c(data.long$polyps00, rep(0,n5+n2)), 
                               l.Adenoma = c(data.long$adenoma00, rep(0,n5+n2)), 
                               z.Time = c(rep(0,n4), log(data.long.z$gap), rep(0,n2)), # related to log(mu): negative binomial component for longitudinal outcomes
                               intercept.z = c(rep(0,n4), rep(1,n5), rep(0,n2)),                               
                               z.Sex = c(rep(0,n4), data.long.z$sex01, rep(0,n2)),
                               z.Age = c(rep(0,n4), data.long.z$age00_trans, rep(0,n2)),
                               z.Pcrc1age = c(rep(0,n4), data.long.z$pcrc1age_trans, rep(0,n2)),
                               z.Cumgap = c(rep(0,n4), data.long.z$cumgap, rep(0,n2)),
                               z.Polyps = c(rep(0,n4), data.long.z$polyps00, rep(0,n2)),
                               z.Adenoma = c(rep(0,n4), data.long.z$adenoma00, rep(0,n2)),
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
  f(d.YL, model="iid", n=n2, hyper=list(prec = list(param = c(5, 0.5)))) +
  f(l.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.YL, copy="d.YL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(d.ZL, model="iid", n=n3, hyper=list(prec = list(param = c(5, 0.5)))) +
  f(l.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(z.ZL, copy="d.ZL", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) 

inla.model1 <- inla(formula1, family = c("binomial","zeroinflatednbinomial0","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE,config = TRUE),
                    control.family = list(list(),list(hyper = list(prob = list(initial = -10,fixed = TRUE))),list(variant=1)),
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
  Pre = 1.67, Running = 0.431, Post = 0.0265, Total = 2.12 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant mode kld
mu          -2.852 0.551     -3.931   -2.852     -1.773   NA   0
intercept.l  1.158 0.391      0.391    1.158      1.925   NA   0
l.Sex       -0.024 0.380     -0.769   -0.024      0.721   NA   0
l.Age       -0.126 0.167     -0.454   -0.126      0.201   NA   0
l.Pcrc1age  -0.122 0.549     -1.197   -0.122      0.954   NA   0
l.Cumgap    -0.117 0.041     -0.197   -0.117     -0.037   NA   0
intercept.z -1.088 0.553     -2.173   -1.088     -0.004   NA   0
z.Sex        1.014 0.512      0.010    1.014      2.018   NA   0
z.Age        0.474 0.242     -0.001    0.474      0.949   NA   0
z.Pcrc1age  -0.790 0.585     -1.937   -0.790      0.357   NA   0
z.Cumgap    -0.297 0.080     -0.454   -0.297     -0.139   NA   0
d.Sex        2.719 0.646      1.453    2.719      3.985   NA   0
d.Adenoma   -0.798 0.843     -2.449   -0.798      0.854   NA   0
d.Age        1.182 0.257      0.679    1.182      1.686   NA   0
d.Pcrc1age  -0.260 0.725     -1.681   -0.260      1.161   NA   0

Random effects:
  Name	  Model
d.YL IID model
d.ZL IID model
l.YL Copy
z.YL Copy
l.ZL Copy
z.ZL Copy

Model hyperparameters:
  mean      sd 0.025quant 0.5quant 0.975quant mode
size for nbinomial zero-inflated observations[2] 36.571 106.103      0.992   12.464    222.126   NA
alpha parameter for weibullsurv[3]                0.482   0.050      0.391    0.479      0.589   NA
Precision for d.YL                                4.650   2.908      1.200    3.961     12.196   NA
Precision for d.ZL                                7.204   3.257      2.825    6.557     15.369   NA
Beta for l.YL                                    -1.377   0.583     -2.527   -1.377     -0.229   NA
Beta for z.YL                                     1.673   0.530      0.607    1.681      2.697   NA
Beta for l.ZL                                     1.991   0.538      0.943    1.986      3.062   NA
Beta for z.ZL                                     1.107   0.620     -0.081    1.094      2.359   NA

Deviance Information Criterion (DIC) ...............: 489.13
Deviance Information Criterion (DIC, saturated) ....: -2446.21
Effective number of parameters .....................: 46.25

Watanabe-Akaike information criterion (WAIC) ...: 498.60
Effective number of parameters .................: 49.38

Marginal log-Likelihood:  -329.07 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model1$marginals.fixed)
#[1] "mu"          "intercept.l" "l.Sex"       "l.Age"       "l.Pcrc1age"  "l.Cumgap"    "intercept.z" "z.Sex"       "z.Age"       "z.Pcrc1age" 
[11] "z.Cumgap"    "d.Sex"       "d.Adenoma"   "d.Age"       "d.Pcrc1age" 
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[1]]))
#Mean            20.1181 
Stdev           11.7563 
Quantile  0.025 5.89104 
Quantile  0.25  11.933 
Quantile  0.5   17.302 
Quantile  0.75  25.0815 
Quantile  0.975 50.7141 

names(inla.model1$internal.marginals.hyperpar) 
#[1] "log size for nbinomial zero-inflated observations[2]" "alpha_intern for weibullsurv[3]"                     
[3] "Log precision for d.YL"                               "Log precision for d.ZL"                              
[5] "Beta_intern for l.YL"                                 "Beta_intern for z.YL"                                
[7] "Beta_intern for l.ZL"                                 "Beta_intern for z.ZL"     

# for dispersion
inla.zmarginal(inla.tmarginal(function(x) exp(x),inla.model1$internal.marginals.hyperpar[[1]]))
#Mean            35.1993 
Stdev           79.7807 
Quantile  0.025 0.620167 
Quantile  0.25  4.66697 
Quantile  0.5   12.0325 
Quantile  0.75  32.1258 
Quantile  0.975 217.252 

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[3]]))
#Mean            0.527061 
Stdev           0.158953 
Quantile  0.025 0.287756 
Quantile  0.25  0.412936 
Quantile  0.5   0.501845 
Quantile  0.75  0.613909 
Quantile  0.975 0.907135 

# for family-specific random effects Z
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[4]]))
#Mean            0.39924 
Stdev           0.0859901 
Quantile  0.025 0.255965 
Quantile  0.25  0.337607 
Quantile  0.5   0.3904 
Quantile  0.75  0.451028 
Quantile  0.975 0.592406 

