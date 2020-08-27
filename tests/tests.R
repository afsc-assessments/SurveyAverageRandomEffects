library(testthat)
source('../R/utils.R')

dir.create('re_runs', showWarnings=FALSE)
dir.create('rem_runs', showWarnings=FALSE)

## Recompile models and move to
setwd('../src')
system('admb re -r')
test <- file.copy('re.exe', to='../tests/re_runs/re.exe', overwrite=TRUE)
if(!test) stop("Failed to make new re.exe model")
system('admb rem -r')
test <- file.copy('rem.exe', to='../tests/rem_runs/rem.exe', overwrite=TRUE)
if(!test) stop("Failed to make new rem.exe model")
setwd('../tests')

library(TMB)
dyn.unload(dynlib('../src/re_tmb'))
compile('../src/re_tmb.cpp')
dyn.load(dynlib('../src/re_tmb'))


## RE model tests

test_re_species <- function(d){
  m <- gsub('.dat','',d)
  file.copy(file.path('data', d), to=file.path('re_runs',d), overwrite=TRUE)
  setwd('re_runs')
  on.exit(setwd('..'))
  message("Testing RE model for ", m)
  file.remove('rwout.rep')
  test <- system(paste('re -ind', d), ignore.stdout=TRUE)
  if(!file.exists('rwout.rep')) stop('Failed to run ', m)
  out <- read_re_output('rwout.rep')
  expect_known_output(out, file=paste0('../testthat/_expect_', m))
  ## Get the data for TMB from the RE output
  srv_sd <- out$srv_sd[!is.na(out$srv_sd)]
  srv_est <- out$srv_est[!is.na(out$srv_est)]
  yrs_srv <- out$year[!is.na(out$srv_est)]
  data <- list(yrs=out$year, yrs_srv=yrs_srv,
               yrs_srv_ind=match(yrs_srv, out$year)-1,
               srv_est=srv_est, srv_sd=srv_sd)
  pars <- list(logSdLam=1, biom=out$log_biomass)
  obj <- MakeADFun(data, pars, random='biom')
  obj$env$beSilent()
  opt <- with(obj, nlminb(par, fn, gr))
  adrep <- sdreport(obj)
  results <- cbind(species=m, out, tmb_log_biomass=adrep$value, tmb_SE_log=adrep$sd)
  return(results)
}

spp <- list.files('data', pattern='.dat')
results <- lapply(spp, test_re_species) %>% do.call(rbind,.)

library(dplyr)
library(tidyr)
library(ggplot2)

results <- results %>%
  mutate(log_biomass_relerror=(log_biomass-tmb_log_biomass)/log_biomass,
         SE_log_relerror=(SE_log-tmb_SE_log)/SE_log) %>%
  pivot_longer(cols=c('SE_log_relerror', 'log_biomass_relerror'),
  values_to='relative_error')

ggplot(results, aes(year, relative_error, color=name)) + geom_point()
