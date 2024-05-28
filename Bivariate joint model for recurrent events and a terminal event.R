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

y.recu <- inla.surv(time = c(time.R/max(time.R), rep(NA, n2)), event = c(delta.R.fij, rep(NA, n2)))
y.term <- inla.surv(time = c(rep(NA, n1), time.D/max(time.D)), event = c(rep(NA, n1), delta.D.fi))
y.joint <- list(y.recu, y.term)

linear.covariate <- data.frame(mu = as.factor(c(rep(1,n1),rep(2,n2))),  # exp(-mu1) and exp(-mu2) = scale parameters in weibull distribution for recurrent events and terminal event respectively
                               r.Sex = c(data$sex01, rep(0,n2)), # relative to hazard function for recurrent events
                               r.Adenoma = c(data$adenoma0, rep(0,n2)),
                               r.Age = c(data$age0_trans, rep(0,n2)),
                               r.Pcrc1age = c(data$pcrc1age_trans, rep(0,n2)),
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
  f(d.YR, model="iid", n=n2, hyper=list(prec = list(param = c(2.5, 0.5)))) +
  f(r.YR, copy="d.YR", hyper = list(beta = list(fixed = FALSE, param = c(0,1)))) +
  f(r.ZR, model="iid", n=n3, hyper=list(prec = list(param = c(2.5, 0.5)))) +
  f(d.ZR, copy="r.ZR", hyper = list(beta = list(fixed = FALSE, param = c(0,1))))

inla.model1 <- inla(formula1, family = c("weibullsurv","weibullsurv"),
                    data = data.f, control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE),
                    control.family = list(list(variant=1),list(variant=1)),
                    control.inla = list(int.strategy="eb", cmin = 0), safe = TRUE,control.fixed = list(prec.intercept = 0.1))

### Final Results
summary(inla.model1)
#Time used:
Pre = 1.24, Running = 1.26, Post = 0.0312, Total = 2.53 
Fixed effects:
  mean    sd 0.025quant 0.5quant 0.975quant   mode kld
mu1         1.424 0.148      1.133    1.424      1.715  1.424   0
mu2        -3.289 0.541     -4.349   -3.289     -2.229 -3.289   0
r.Sex       0.282 0.150     -0.012    0.282      0.575  0.282   0
r.Adenoma   0.176 0.156     -0.129    0.176      0.481  0.176   0
r.Age       0.013 0.061     -0.108    0.013      0.133  0.013   0
r.Pcrc1age -0.305 0.261     -0.815   -0.305      0.206 -0.305   0
d.Sex       3.376 0.699      2.007    3.376      4.746  3.376   0
d.Adenoma  -1.944 0.911     -3.730   -1.944     -0.158 -1.944   0
d.Age       1.244 0.295      0.665    1.244      1.822  1.244   0
d.Pcrc1age -0.976 0.795     -2.534   -0.976      0.582 -0.976   0

Random effects:
  Name	  Model
d.YR IID model
r.ZR IID model
r.YR Copy
d.ZR Copy

Model hyperparameters:
  mean    sd 0.025quant 0.5quant 0.975quant  mode
alpha parameter for weibullsurv     1.201 0.064      1.078    1.199      1.331 1.197
alpha parameter for weibullsurv[2]  0.617 0.022      0.577    0.616      0.663 0.613
Precision for d.YR                  0.157 0.054      0.075    0.150      0.284 0.135
Precision for r.ZR                  5.975 1.839      3.165    5.708     10.336 5.209
Beta for r.YR                       0.206 0.049      0.110    0.206      0.301 0.207
Beta for d.ZR                      -0.016 0.622     -1.272   -0.006      1.177 0.039

Deviance Information Criterion (DIC) ...............: -319.67
Deviance Information Criterion (DIC, saturated) ....: 806.33
Effective number of parameters .....................: 100.30

Watanabe-Akaike information criterion (WAIC) ...: -268.51
Effective number of parameters .................: 120.27

Marginal log-Likelihood:  79.02 
CPO, PIT is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

# scale parameter in Weibull baseline hazard function for recurrent event and a terminal event
# b_R = exp(-intercept_R) and b_D = exp(-intercept_D)
names(inla.model1$marginals.fixed)
#[1] "mu1"        "mu2"        "r.Sex"      "r.Adenoma"  "r.Age"      "r.Pcrc1age" "d.Sex"      "d.Adenoma"  "d.Age"      "d.Pcrc1age"
inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[1]]))
#Mean            0.243388 
Stdev           0.0360705 
Quantile  0.025 0.180158 
Quantile  0.25  0.217792 
Quantile  0.5   0.240696 
Quantile  0.75  0.266008 
Quantile  0.975 0.321573 

inla.zmarginal(inla.tmarginal(function(x) exp(-x),inla.model1$marginals.fixed[[2]]))
#Mean            30.9667 
Stdev           17.7406 
Quantile  0.025 9.28893 
Quantile  0.25  18.5838 
Quantile  0.5   26.7708 
Quantile  0.75  38.5569 
Quantile  0.975 77.0095 

names(inla.model1$internal.marginals.hyperpar) 
#[1] "alpha_intern for weibullsurv"    "alpha_intern for weibullsurv[2]" "Log precision for d.YR"          "Log precision for r.ZR"         
[5] "Beta_intern for r.YR"            "Beta_intern for d.ZR"      

# for standard deviation and correlation coefficients of random effects
# for subject-specific random effects Y
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[3]]))
#Mean            2.63222 
Stdev           0.451492 
Quantile  0.025 1.88125 
Quantile  0.25  2.30944 
Quantile  0.5   2.58282 
Quantile  0.75  2.90371 
Quantile  0.975 3.64923 

# for family-specific random effects Z
inla.zmarginal(inla.tmarginal(function(x) sqrt(1/exp(x)),inla.model1$internal.marginals.hyperpar[[4]]))
#Mean            0.423129 
Stdev           0.0634436 
Quantile  0.025 0.311782 
Quantile  0.25  0.378108 
Quantile  0.5   0.418464 
Quantile  0.75  0.462954 
Quantile  0.975 0.5605 

inla.model1$dic$dic
#[1]  -319.6706
inla.model1$dic$family.dic
#[1]  -287.4758  -32.1948
