## Quick file to explore some things we can do with the TMB
## version of the model.

library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(TMB)
library(tmbstan)
## dyn.unload(dynlib('../src/re_tmb'))
## file.remove('../src/re_tmb.dll')
compile('src/re_tmb.cpp')
dyn.load(dynlib('src/re_tmb'))
source('R/utils.R')

out <- read_re_output('tests/re_runs/rwout.rep')
srv_sd <- out$srv_sd[!is.na(out$srv_sd)]
srv_est <- out$srv_est[!is.na(out$srv_est)]
yrs_srv <- out$year[!is.na(out$srv_est)]
## Convert SD back to CV which then gets converted internally
## back to SD
srv_cv <- sqrt(exp(srv_sd^2)-1)
data <- list(yrs=out$year, yrs_srv=yrs_srv,
             yrs_srv_ind=match(yrs_srv, out$year)-1,
             srv_est=srv_est, srv_cv=srv_cv)
pars <- list(logSdLam=1, biom=out$log_biomass)
obj <- MakeADFun(data, pars, random='biom')
obj$env$beSilent()
obj$fn()
opt <- with(obj, nlminb(par, fn, gr))
adrep <- sdreport(obj)

## Get MLE estimates
biom.mle <- with(adrep,
                 data.frame(year=data$yrs,
                            med=value,
                            upr=value+1.96*sd,
                            lwr=value-1.96*sd,
                            type='mle'))

fit <- tmbstan(obj)
post <- as.data.frame(fit)

## Compare MCMC vs MLE of biomass
biom.mcmc <- list()
for(i in 1:nrow(fit)){
  biom.mcmc[[i]] <- data.frame(year=data$yrs,
                               biom=obj$report(par=post[i,-ncol(post)])$biom)
}
biom.mcmc <- bind_rows(biom.mcmc) %>%
  group_by(year) %>%
  summarize(lwr=quantile(biom,.025),
            upr=quantile(biom, .975), med=median(biom)) %>%
  mutate(type='mcmc')
biom <- bind_rows(biom.mcmc, biom.mle)
g <- ggplot(biom, aes(year, med, ymin=lwr, ymax=upr, color=type, fill=type)) +
  geom_ribbon(alpha=.5) + geom_line() + labs(y='log biomass')
ggsave('mcmc_vs_mle.png', g, width=7, height=5)

## And the estimate of logSdLam
hist(post[,1], freq=FALSE)
x <- seq(min(post[,1]), max(post[1,]), len=500)
lines(x, y=dnorm(x, adrep$par.fixed[1], sqrt(adrep$cov.fixed[1,1])), lwd=2)

## Check the LA
fit1 <- tmbstan(obj, chains=4, laplace=FALSE, iter=10000,
                thin=10)
fit2 <- tmbstan(obj, chains=4, laplace=TRUE, iter=10000,
                thin=10)
post1 <- as.data.frame(fit1)[,1]
post2 <- as.data.frame(fit2)[,1]
monitor(post1, print=FALSE)[,'n_eff']
monitor(post2, print=FALSE)[,'n_eff']

png('qq_la_check.png', width=4, height=4, units='in', res=300)
qqplot(post1,post2); abline(0,1)
dev.off()
