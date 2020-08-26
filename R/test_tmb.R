
library(TMB)

dyn.unload(dynlib('../src/re_tmb'))
compile('../src/re_tmb.cpp')
dyn.load(dynlib('../src/re_tmb'))


## Quick example from Bogoslof pollock
styr <- 1990
endyr <- 2013
yrs_srv <- c(1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2005,2006,2007,2009,2012)
srv_est <- c(1289006,940198,635405,490077,1104124,682277,392402,492396,475311,301402,232170,225712,197851,253459,240000,292000,110000,67500)
srv_se <- c(154680.72, 188039.6, 57186.45, 58809.24, 121453.64, 136455.4, 54936.28, 93555.24, 104568.42, 42196.28, 23217, 27085.44, 43527.22, 43088.03, 28800, 35040, 20900, 12825)
srv_cv  <- sqrt(log(1+(srv_se/srv_est)^2))
## Convert to TMB inputs
data <- list(styr=styr, endyr=endyr,
             yrs_srv=yrs_srv,
             yrs_srv_ind=match(yrs_srv, styr:endyr)-1,
             srv_est=srv_est,
             srv_cv=srv_cv)
pars <- list(logSdLam=1, biom=rnorm(n=endyr-styr+1))

obj <- MakeADFun(data, pars, random='biom')
obj$fn()

opt <- with(obj, nlminb(par, fn, gr))
rep <- obj$report()
adrep <- sdreport(obj)
tmb <- data.frame(logbiomass=adrep$value, logbiomass.sd=adrep$sd)


## Compare to RE model in ADMB (copied manually from output file)
admb <- data.frame(biomsd=c(14.0357, 14.0357, 13.7351, 13.3735, 13.1927,
                    13.8178, 13.417, 12.9634, 13.0577, 12.9746,
                    12.6295, 12.3717, 12.3243, 12.2528, 12.3367,
                    12.4205, 12.4037, 12.5237, 12.101, 11.6784,
                    11.508, 11.3377, 11.1673, 11.1673),
           biomsd.sd=c(0.375136, 0.115255, 0.158327, 0.085484,
                       0.112331, 0.106645, 0.158958, 0.127109,
                       0.15629, 0.172169, 0.124107, 0.0936449,
                       0.109628, 0.177765, 0.279996, 0.148073,
                       0.109438, 0.113246, 0.273526, 0.172203,
                       0.321117, 0.322879, 0.181822, 0.400628))
round((admb$biomsd-tmb$logbiomass)/admb$biomsd,4)
round((admb$biomsd.sd-tmb$logbiomass.sd)/admb$biomsd.sd,4)

plot(x=styr:endyr, y=rep$biom, type='b')
points(data$yrs_srv, log(data$srv_est), col=2)
points(styr:endyr, admb$biomsd, col=3,cex=.5)


## Check LA approximation with tmbstan
library(tmbstan)
mcmcfit1 <- tmbstan(obj, laplace=TRUE)
obj2 <- MakeADFun(data, pars)
mcmcfit2 <- tmbstan(obj2)

post1 <- as.data.frame(mcmcfit1)
post2 <- as.data.frame(mcmcfit2)
qqplot(post1[,1], post2[,1]); abline(0,1, col=2)
