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

## RE model tests
file.copy('data/AI_Kamchatka.dat', to='re_runs/AI_Kamchatka.dat', overwrite=TRUE)
setwd('re_runs')
system('re -ind AI_Kamchatka.dat')
out <- read_re_output('rwout.rep')
setwd('..')

expect_known_output(out, file='testthat/_expect_AI_Kamchatka')
