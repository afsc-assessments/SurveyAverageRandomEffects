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
### this and checks it is the same.

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
results <- test_re_species("bog.dat")
results <- results %>%
  mutate(log_biomass_relerror=(log_biomass-tmb_log_biomass)/log_biomass,
         SE_log_relerror=(SE_log-tmb_SE_log)/SE_log) %>%
  pivot_longer(cols=c('SE_log_relerror', 'log_biomass_relerror'),
  values_to='relative_error')
ggplot(results, aes(year, relative_error, color=name)) + geom_point()
