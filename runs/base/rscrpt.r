sims_pop_tru=read.table("poptrusrv.out",header=T)
sims_wp_tru=read.table("wptrusrv.out",header=T)
ln_sims_wp=read.table("ln_sims_wp_2013.rep")
sims_pop=read.table("sims_pop_2013.rep")
ln_sims_pop=read.table("ln_sims_pop_2013.rep")
sims_wp=read.table("sims_wp_2013.rep")
scen=read.table("scen.dat",header=F)
# Relative biomass for ln pollock
ln_rb_wp=ln_sims_wp[,55]/sims_wp_tru$X54-1
# Relative biomass for normal pollock
rb_wp=sims_wp[,55]/sims_wp_tru$X54-1
# Relative biomass for log-normal pop
ln_rb_pop=ln_sims_pop[,55]/sims_pop_tru$X54-1
# Relative biomass for normal pop
rb_pop=sims_pop[,55]/sims_pop_tru$X54-1

# summarize them
mn_lrb_wp=tapply(ln_rb_wp,scen$V1,FUN=mean)
mn_rb_wp=tapply(rb_wp,scen$V1,FUN=mean)

# Plot all of them
# lognormal
boxplot(ln_rb_wp~scen$V1,ylim=c(-.5,.5),col="salmon",
        xaxt="n",main="Pollock, log-Normal")
abline(h=0)
axis(1,at=1:36,rownames(mn_lrb_wp),las=2,cex.axis=.8)

# normal
boxplot(rb_wp~scen$V1,ylim=c(-.5,.5),col="salmon",
        xaxt="n",main="Pollock")
axis(1,at=1:36,rownames(mn_lrb_wp),las=2,cex.axis=.8)
abline(h=0)


plot(mn_rb_wp/mn_lrb_wp,pch=19,ylim=c(-2,2))
abline(h=0)

par(mai=c(2,1.2,1,.5))
plot(mn_lrb_wp,ylim=c(-.3,.3),pch=19,xaxt="n",xlab="",ylab="Relative Biomass Error")
points(mn_rb_wp,pch=1)
abline(h=0)
axis(1,at=1:36,rownames(x),las=2,cex.axis=.8)
legend(1,.25,legend=c("log-normal","normal"),pch=c(19,1))

boxplot(ln_rb_wp[1:1200],rb_wp[1:1200],xaxt="n",main="")
boxplot(ln_rb_wp[1:1200],ln_rb_wp[2401:3600],rb_wp[1:1200],rb_wp[2401:3600],
        ylim=c(-.7,.7),xaxt="n",main="Pollock-like",col="gold")
abline(h=0,lwd=3,lty=2,col="blue")
axis(1,at=1:4,c("DownUpF LNormal","UpDownF LNormal","DownUpF Normal","UpDownF Normal"))

