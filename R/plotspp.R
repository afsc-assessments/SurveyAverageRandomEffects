R
source("~/OneDrive/Models/Atka/2018/R/prelims.R")
source("read-admb.R")
library(ggplot2)
library(ggthemes)
library(tidyr)
library(dplyr)

b <- read_rep("arc/bog.rep")
t <- read_rep("arc/atka_t.rep")
e <- read_rep("arc/atka_e.rep")
w <- read_rep("arc/atka_w.rep")
c <- read_rep("arc/atka_c.rep")
M <- b
plot_spp <- function(M,main="shit",ymax=800){
  t.df <- data.frame(Year=M$yrs, lci=M$LCI/1e3, Biomass=M$biomA/1e3, uci=M$UCI/1e3, l90=M$low90th/1e3,u90=M$upp90th/1e3)
  d.df <- data.frame(Year=M$yrs_srv, lci=(M$srv_est-M$srv_est*M$srv_sd)/1e3, uci=(M$srv_est+M$srv_est*M$srv_sd)/1e3, Biomass=M$srv_est/1e3)
  p <- ggplot(t.df, aes(x=Year,y=Biomass)) + geom_ribbon(aes(ymin=l90,ymax=u90),color="brown",fill="salmon",alpha=.4) + 
    geom_line() + theme_few() + ggtitle(main) + coord_cartesian(ylim=c(0, ymax)) + ylab("Biomass (kt)") +  
    geom_point(data=d.df,aes(x=Year,y=Biomass),size=2,color="blue") + geom_errorbar(data=d.df,aes(x=Year,ymin=lci,ymax=uci),color="blue",width=.3) 
    print(paste0("Random effects biomass for ",t.df[dim(t.df)[1],1]," in ",main," is ",round(t.df[dim(t.df)[1],3],0), " thousand t") )
    return(p)
}
ymax=1500
    #geom_ribbon(aes(ymin=lci,ymax=uci),fill="salmon",alpha=.4) + 
plot_spp(b,"Bogoslof Island",ymax=1800)
plot_spp(e,"Eastern Aleutians")
plot_spp(c,"Central Aleutians")
plot_spp(w,"Western Aleutians")
plot_spp(w,"Western Aleutians")

library(scales)
show_col(fivethirtyeight_pal)
p <- ggplot(df, aes(trt, resp, colour = group))
     p + geom_linerange(aes(ymin = lower, ymax = upper))
     p + geom_pointrange(aes(ymin = lower, ymax = upper))
     p + geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.2)
     p + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
     
     # Draw lines connecting group means
     p +
       geom_line(aes(group = group)) +
       geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
     
     # If you want to dodge bars and errorbars, you need to manually
     # specify the dodge width
     p <- ggplot(df, aes(trt, resp, fill = group))
     p +
      geom_col(position = "dodge") +
      geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.25)
     
     # Because the bars and errorbars have different widths
     # we need to specify how wide the objects we are dodging are
     dodge <- position_dodge(width=0.9)
     p +
       geom_col(position = dodge) +
       geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = 0.25)
     
     # When using geom_errorbar() with position_dodge2(), extra padding will be
     # needed between the error bars to keep them aligned with the bars.
     p +
     geom_col(position = "dodge2") +
     geom_errorbar(
       aes(ymin = lower, ymax = upper),
       position = position_dodge2(width = 0.5, padding = 0.5)
