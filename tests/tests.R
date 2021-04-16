### This file tests a set of model inputs against previous runs
### to check that the results are identical. This is useful to
### check if either the RE model changes, or if ADMB is
### updated. It is setup to do the univariate RE model but
### not currently the multivariate REM model. The previous
### (expected) results are stored in the folder tests/testthat as
### _expect* files. If there is a need to update these files,
### just delete them and rerun this code and it will recreate
### them. See help for ?testthat for more information.

### In addition I created a TMB model version and this also runs
### this and checks it is the same. There are some minor
### differences so a tolerance of 1e-4 is used. I believe this is
### due to differences in the inner optimizer (i.e., the RE don't
### match exactly given the same fixed effects), but was unable
### to confirm this (see below).

### To run, open this script and set the working directory to the
### 'test' folder then run below.

### Cole Monnahan | April 2021


### Setup testing environment
library(testthat)
source('../R/utils.R')
dir.create('re_runs', showWarnings=FALSE)
dir.create('rem_runs', showWarnings=FALSE)

### Recompile ADMD-RE models and move to testing folders.
setwd('../src')
system('admb -r re ')
test <- file.copy('re.exe', to='../tests/re_runs/re.exe', overwrite=TRUE)
## Non-Windows version
## test <- file.copy('re', to='../tests/re_runs/re', overwrite=TRUE)
if(!test) stop("Failed to make new re.exe model")
system('admb -r rem ')
test <- file.copy('rem.exe', to='../tests/rem_runs/rem.exe', overwrite=TRUE)
##test <- file.copy('rem', to='../tests/rem_runs/rem', overwrite=TRUE)
if(!test) stop("Failed to make new rem.exe model")
setwd('../tests')

### Recompile TMB models
library(TMB)
## dyn.unload(dynlib('../src/re_tmb'))
## file.remove('../src/re_tmb.dll')
compile('../src/re_tmb.cpp')
dyn.load(dynlib('../src/re_tmb'))

### Run tests across all species that have no zeroes. See
### "broken" for those that do and thus fail here.
spp <- list.files('data', pattern='.dat')
results <- lapply(spp, test_re_species) %>% do.call(rbind,.)

### Informal test of TMB vs ADMB
library(dplyr)
library(tidyr)
library(ggplot2)
g <- results %>%
  mutate(log_biomass_re=(log_biomass-tmb_log_biomass)/log_biomass,
         SE_log_re=(SE_log-tmb_SE_log)/SE_log) %>%
  pivot_longer(cols=c('SE_log_re', 'log_biomass_re'),
               values_to='relative_error') %>%
  ggplot(aes(species, abs(relative_error), col=log_biomass)) + geom_point() +
  facet_wrap('name', ncol=1) + scale_y_log10() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x=NA, y='Absolute relative error (log space)',
       title='ADMB-RE vs TMB Random Effects Model Differences')
ggsave('TMB_vs_ADMB_comparison.png', g, width=8, height=6)


## ## There are some small differences. I think this is due to the
## ## inner optimizer being different. Below is some old code that
## ## was trying to compare the two. It doesn't work but we
## ## might want to revisit it later.
## compile('../src/re_tmb.cpp')
## dyn.load(dynlib('../src/re_tmb'))
## set.seed(3234)
## yrs <- 1:30
## est <- sort(rnorm(30, 1e12, 1e12*.1))
## ## est <- est + runif(30)
## CV <- rep(.1, 30)
## plot(est)
## ## Fit in ADMB first to get MLE
## write_re_input('re_runs/test.dat', years=yrs, biomass=est, CV=CV,
##                tag='test to compare ADMB and TMB')
## setwd('re_runs/')
## system('re -ind test.dat')
## out <- read_re_output()
## setwd('..')
## par.mle <- as.numeric(readLines('re_runs/re.par')[3])
## ## Use MLE as starting value for TMB
## data <- list(yrs=yrs, yrs_srv=yrs,
##              yrs_srv_ind=yrs-1,
##              srv_est=est, srv_cv=CV)
## pars <- list(logSdLam=par.mle, biom=log(est))
## obj <- MakeADFun(data, pars, random='biom')
## obj$env$beSilent()
## ## opt <- with(obj, nlminb(par, fn, gr))
## adrep <- sdreport(obj)
## results <- test_re_species('test.dat')
## results %>%
##   transmute(RE_est=(tmb_log_biomass-log_biomass)/log_biomass,
##             RE_sd=(tmb_SE_log-SE_log)/SE_log) %>% head(3)
