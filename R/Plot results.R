#============================================================================================
#============= Get Working Directories
#============================================================================================

path<-getwd()

pathD<-paste(path,"/Plot data",sep="")

#============================================================================================
#============= Read in data/set up sensitivity analysis parameters,results matrices
#============================================================================================

CGOA_re<-read.csv(paste(pathD,"/CGOA_re.csv",sep=""))
CGOA_re_LL<-read.csv(paste(pathD,"/CGOA_re_LL.csv",sep=""))
CGOA_re_17<-read.csv(paste(pathD,"/CGOA_re_17.csv",sep=""))
CGOA_srv<-read.csv(paste(pathD,"/CGOA_srv.csv",sep=""))
CGOA_srv_LL<-read.csv(paste(pathD,"/CGOA_srv_LL.csv",sep=""))

WGOA_re<-read.csv(paste(pathD,"/WGOA_re.csv",sep=""))
WGOA_re_LL<-read.csv(paste(pathD,"/WGOA_re_LL.csv",sep=""))
WGOA_re_17<-read.csv(paste(pathD,"/WGOA_re_17.csv",sep=""))
WGOA_srv<-read.csv(paste(pathD,"/WGOA_srv.csv",sep=""))
WGOA_srv_LL<-read.csv(paste(pathD,"/WGOA_srv_LL.csv",sep=""))

EGOA_re<-read.csv(paste(pathD,"/EGOA_re.csv",sep=""))
EGOA_re_LL<-read.csv(paste(pathD,"/EGOA_re_LL.csv",sep=""))
EGOA_re_17<-read.csv(paste(pathD,"/EGOA_re_17.csv",sep=""))
EGOA_srv<-read.csv(paste(pathD,"/EGOA_srv.csv",sep=""))
EGOA_srv_LL<-read.csv(paste(pathD,"/EGOA_srv_LL.csv",sep=""))

GOA_re<-read.csv(paste(pathD,"/GOA_re.csv",sep=""))
GOA_re_LL<-read.csv(paste(pathD,"/GOA_re_LL.csv",sep=""))
GOA_re_17<-read.csv(paste(pathD,"/GOA_re_17.csv",sep=""))
GOA_srv<-read.csv(paste(pathD,"/GOA_srv.csv",sep=""))
GOA_srv_LL<-read.csv(paste(pathD,"/GOA_srv_LL.csv",sep=""))

#============================================================================================
#============= Plot Total trawl survey (for SAFE)
#============================================================================================

par(mar=c(3,6,0.1,1))

x<-plot(GOA_re$yrs,GOA_re$biomA,type="l",lwd=0,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,120000),cex.axis=1.25,col="white")
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est+1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est-1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
lines(GOA_srv$yrs_srv,GOA_srv$srv_est,lty=3)
pch_plot<-c(16,16,16,16,16,16,21,16,16,16,16,16,16,16,16)
points(GOA_srv$yrs_srv,GOA_srv$srv_est,pch=pch_plot,col="aquamarine4",bg="white",cex=2)
mtext("Biomass (mt)",side=2,line=4.5,cex=1.25)
axis(1,cex.axis=1.5)


#============================================================================================
#============= Plot Total fit (bottom trawl survey with LL survey)
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,0),4,1,byrow = TRUE),heights=c(0.1,1,1,0.2))
par(mar=c(0.1,9,0.1,1))
#options(scipen=999)

x<-plot(GOA_re$yrs,GOA_re$biomA,type="l",lwd=3,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,120000),cex.axis=1.5)
polygon(c(GOA_re$yrs,sort(GOA_re$yrs,decreasing=TRUE)),c(GOA_re$UCI,rev(GOA_re$LCI)),col="light grey",border=NA)
lines(GOA_re_17$yrs,GOA_re_17$biomA,lwd=2)
lines(GOA_re$yrs,GOA_re$biomA,lwd=3,col="orangered")
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est+1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est-1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
pch_plot<-c(16,16,16,16,16,16,21,16,16,16,16,16,16,16,16)
points(GOA_srv$yrs_srv,GOA_srv$srv_est,pch=pch_plot,col="aquamarine4",bg="white",cex=2)
legend(1983,120777,legend=c("Bottom trawl survey biomass","2017.0","2019.2b"),col=c("aquamarine4","black","orangered"),lty=c(NA,1,1),lwd=c(NA,2,3),pch=c(16,NA,NA),cex=1.5)
mtext("Biomass (mt)",side=2,line=6,cex=1.25)

x<-plot(GOA_re_LL$yrs,GOA_re_LL$biomA,type="l",lwd=3,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,80000),cex.axis=1.5)
polygon(c(GOA_re_LL$yrs,sort(GOA_re_LL$yrs,decreasing=TRUE)),c(GOA_re_LL$UCI,rev(GOA_re_LL$LCI)),col="light grey",border=NA)
lines(GOA_re_LL$yrs,GOA_re_LL$biomA,lwd=3,col="orangered")
arrows(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est+1.96*GOA_srv_LL$srv_est*GOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est-1.96*GOA_srv_LL$srv_est*GOA_srv_LL$srv_sd,angle=90,length=0.05)
#pch_plot<-c(16,16,21,21,21,16,21,21,16,16,16,21,21,16,21)
points(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,pch=16,col="darkorange4",bg="white",cex=2)
legend(1983,80000,legend=c("Longline survey RPW","2019.2b"),col=c("darkorange4","orangered"),lty=c(NA,1),lwd=c(NA,3),pch=c(16,NA),cex=1.5)
mtext("Relative Population Weight (RPW)",side=2,line=6,cex=1.25)
axis(1,cex.axis=1.5)



#============================================================================================
#============= Plot Total fit (bottom trawl survey with LL survey) for SAFE
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,0),4,1,byrow = TRUE),heights=c(0.1,1,1,0.2))
par(mar=c(0.1,9,0.1,1))
#options(scipen=999)

x<-plot(GOA_re$yrs,GOA_re$biomA,type="l",lwd=3,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,120000),cex.axis=1.5)
polygon(c(GOA_re$yrs,sort(GOA_re$yrs,decreasing=TRUE)),c(GOA_re$UCI,rev(GOA_re$LCI)),col="light grey",border=NA)
lines(GOA_re$yrs,GOA_re$biomA,lwd=3)
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est+1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
arrows(GOA_srv$yrs_srv,GOA_srv$srv_est,GOA_srv$yrs_srv,GOA_srv$srv_est-1.96*GOA_srv$srv_est*GOA_srv$srv_sd,angle=90,length=0.05)
pch_plot<-c(16,16,16,16,16,16,21,16,16,16,16,16,16,16,16)
points(GOA_srv$yrs_srv,GOA_srv$srv_est,pch=pch_plot,col="aquamarine4",bg="white",cex=2)
mtext("Biomass (mt)",side=2,line=6,cex=1.25)

x<-plot(GOA_re_LL$yrs,GOA_re_LL$biomA,type="l",lwd=3,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,80000),cex.axis=1.5)
polygon(c(GOA_re_LL$yrs,sort(GOA_re_LL$yrs,decreasing=TRUE)),c(GOA_re_LL$UCI,rev(GOA_re_LL$LCI)),col="light grey",border=NA)
lines(GOA_re_LL$yrs,GOA_re_LL$biomA,lwd=3)
arrows(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est+1.96*GOA_srv_LL$srv_est*GOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est-1.96*GOA_srv_LL$srv_est*GOA_srv_LL$srv_sd,angle=90,length=0.05)
points(GOA_srv_LL$yrs_srv,GOA_srv_LL$srv_est,pch=16,col="darkorange4",bg="white",cex=2)
mtext("Relative Population Weight (RPW)",side=2,line=6,cex=1.25)
axis(1,cex.axis=1.5)


#============================================================================================
#============= Plot Re fit to bottom trawl by region - same scale
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,3,0),5,1,byrow = TRUE),heights=c(0.1,1,1,1,0.2))
par(mar=c(0.1,7,0.1,1))
#options(scipen=999)

pch_plot<-c(16,16,16,16,16,16,21,16,16,16,16,16,16,16,16)

x<-plot(WGOA_re$yrs,WGOA_re$biomA,type="l",lwd=3,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,20000))
polygon(c(WGOA_re$yrs,sort(WGOA_re$yrs,decreasing=TRUE)),c(WGOA_re$UCI,rev(WGOA_re$LCI)),col="light grey",border=NA)
lines(WGOA_re_17$yrs,WGOA_re_17$biomA,lwd=2)
lines(WGOA_re$yrs,WGOA_re$biomA,lwd=3,col="orangered")
arrows(WGOA_srv$yrs_srv,WGOA_srv$srv_est,WGOA_srv$yrs_srv,WGOA_srv$srv_est+1.96*WGOA_srv$srv_est*WGOA_srv$srv_sd,angle=90,length=0.05)
arrows(WGOA_srv$yrs_srv,WGOA_srv$srv_est,WGOA_srv$yrs_srv,WGOA_srv$srv_est-1.96*WGOA_srv$srv_est*WGOA_srv$srv_sd,angle=90,length=0.05)
points(WGOA_srv$yrs_srv,WGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,17770,"WGOA",cex=1.5)
legend(1987,20000,legend=c("Bottom trawl survey biomass","2017.0","2019.2b"),col=c("aquamarine4","black","orangered"),lty=c(NA,1,1),lwd=c(NA,2,3),pch=c(16,NA,NA),cex=1.5)


x<-plot(CGOA_re$yrs,CGOA_re$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,50000))
polygon(c(CGOA_re$yrs,sort(CGOA_re$yrs,decreasing=TRUE)),c(CGOA_re$UCI,rev(CGOA_re$LCI)),col="light grey",border=NA)
lines(CGOA_re_17$yrs,CGOA_re_17$biomA,lwd=2)
lines(CGOA_re$yrs,CGOA_re$biomA,lwd=3,col="orangered")
arrows(CGOA_srv$yrs_srv,CGOA_srv$srv_est,CGOA_srv$yrs_srv,CGOA_srv$srv_est+1.96*CGOA_srv$srv_est*CGOA_srv$srv_sd,angle=90,length=0.05)
arrows(CGOA_srv$yrs_srv,CGOA_srv$srv_est,CGOA_srv$yrs_srv,CGOA_srv$srv_est-1.96*CGOA_srv$srv_est*CGOA_srv$srv_sd,angle=90,length=0.05)
points(CGOA_srv$yrs_srv,CGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,47770,"CGOA",cex=1.5)
mtext("Biomass (mt)",side=2,line=4.5,cex=1.25)

x<-plot(EGOA_re$yrs,EGOA_re$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,70000))
polygon(c(EGOA_re$yrs,sort(EGOA_re$yrs,decreasing=TRUE)),c(EGOA_re$UCI,rev(EGOA_re$LCI)),col="light grey",border=NA)
lines(EGOA_re_17$yrs,EGOA_re_17$biomA,lwd=2)
lines(EGOA_re$yrs,EGOA_re$biomA,lwd=3,col="orangered")
arrows(EGOA_srv$yrs_srv,EGOA_srv$srv_est,EGOA_srv$yrs_srv,EGOA_srv$srv_est+1.96*EGOA_srv$srv_est*EGOA_srv$srv_sd,angle=90,length=0.05)
arrows(EGOA_srv$yrs_srv,EGOA_srv$srv_est,EGOA_srv$yrs_srv,EGOA_srv$srv_est-1.96*EGOA_srv$srv_est*EGOA_srv$srv_sd,angle=90,length=0.05)
points(EGOA_srv$yrs_srv,EGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,67770,"EGOA",cex=1.5)

axis(1,cex.axis=1.25)



#============================================================================================
#============= Plot Re fit to bottom trawl by region - same scale (for SAFE)
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,3,0),5,1,byrow = TRUE),heights=c(0.1,1,1,1,0.2))
par(mar=c(0.1,7,0.1,1))
#options(scipen=999)

pch_plot<-c(16,16,16,16,16,16,21,16,16,16,16,16,16,16,16)

x<-plot(WGOA_re$yrs,WGOA_re$biomA,type="l",lwd=3,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,30000))
polygon(c(WGOA_re$yrs,sort(WGOA_re$yrs,decreasing=TRUE)),c(WGOA_re$UCI,rev(WGOA_re$LCI)),col="light grey",border=NA)
lines(WGOA_re$yrs,WGOA_re$biomA,lwd=3)
arrows(WGOA_srv$yrs_srv,WGOA_srv$srv_est,WGOA_srv$yrs_srv,WGOA_srv$srv_est+1.96*WGOA_srv$srv_est*WGOA_srv$srv_sd,angle=90,length=0.05)
arrows(WGOA_srv$yrs_srv,WGOA_srv$srv_est,WGOA_srv$yrs_srv,WGOA_srv$srv_est-1.96*WGOA_srv$srv_est*WGOA_srv$srv_sd,angle=90,length=0.05)
points(WGOA_srv$yrs_srv,WGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,27770,"WGOA",cex=1.5)


x<-plot(CGOA_re$yrs,CGOA_re$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,70000))
polygon(c(CGOA_re$yrs,sort(CGOA_re$yrs,decreasing=TRUE)),c(CGOA_re$UCI,rev(CGOA_re$LCI)),col="light grey",border=NA)
lines(CGOA_re$yrs,CGOA_re$biomA,lwd=3)
arrows(CGOA_srv$yrs_srv,CGOA_srv$srv_est,CGOA_srv$yrs_srv,CGOA_srv$srv_est+1.96*CGOA_srv$srv_est*CGOA_srv$srv_sd,angle=90,length=0.05)
arrows(CGOA_srv$yrs_srv,CGOA_srv$srv_est,CGOA_srv$yrs_srv,CGOA_srv$srv_est-1.96*CGOA_srv$srv_est*CGOA_srv$srv_sd,angle=90,length=0.05)
points(CGOA_srv$yrs_srv,CGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,67770,"CGOA",cex=1.5)
mtext("Biomass (mt)",side=2,line=4.5,cex=1.25)

x<-plot(EGOA_re$yrs,EGOA_re$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,70000))
polygon(c(EGOA_re$yrs,sort(EGOA_re$yrs,decreasing=TRUE)),c(EGOA_re$UCI,rev(EGOA_re$LCI)),col="light grey",border=NA)
lines(EGOA_re$yrs,EGOA_re$biomA,lwd=3)
arrows(EGOA_srv$yrs_srv,EGOA_srv$srv_est,EGOA_srv$yrs_srv,EGOA_srv$srv_est+1.96*EGOA_srv$srv_est*EGOA_srv$srv_sd,angle=90,length=0.05)
arrows(EGOA_srv$yrs_srv,EGOA_srv$srv_est,EGOA_srv$yrs_srv,EGOA_srv$srv_est-1.96*EGOA_srv$srv_est*EGOA_srv$srv_sd,angle=90,length=0.05)
points(EGOA_srv$yrs_srv,EGOA_srv$srv_est,pch=pch_plot,cex=1.75,col="aquamarine4",bg="white")
text(1985,67770,"EGOA",cex=1.5)

axis(1,cex.axis=1.25)



#============================================================================================
#============= Plot Re fit to LL by region - same scale
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,3,0),5,1,byrow = TRUE),heights=c(0.1,1,1,1,0.2))
par(mar=c(0.1,7,0.1,1))

x<-plot(WGOA_re_LL$yrs,WGOA_re_LL$biomA,type="l",lwd=3,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,20000))
polygon(c(WGOA_re_LL$yrs,sort(WGOA_re_LL$yrs,decreasing=TRUE)),c(WGOA_re_LL$UCI,rev(WGOA_re_LL$LCI)),col="light grey",border=NA)
lines(WGOA_re_LL$yrs,WGOA_re_LL$biomA,lwd=3,col="orangered")
arrows(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est+1.96*WGOA_srv_LL$srv_est*WGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est-1.96*WGOA_srv_LL$srv_est*WGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,17770,"WGOA",cex=1.5)
legend(1984,15000,legend=c("Longline survey RPW","2019.2b"),col=c("darkorange4","orangered"),lty=c(NA,1),lwd=c(NA,3),pch=c(16,NA),cex=1.5)


x<-plot(CGOA_re_LL$yrs,CGOA_re_LL$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,20000))
polygon(c(CGOA_re_LL$yrs,sort(CGOA_re_LL$yrs,decreasing=TRUE)),c(CGOA_re_LL$UCI,rev(CGOA_re_LL$LCI)),col="light grey",border=NA)
lines(CGOA_re_LL$yrs,CGOA_re_LL$biomA,lwd=3,col="orangered")
arrows(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est+1.96*CGOA_srv_LL$srv_est*CGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est-1.96*CGOA_srv_LL$srv_est*CGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,17770,"CGOA",cex=1.5)
mtext("Relative Population Weight (RPW)",side=2,line=4.5,cex=1.25)

x<-plot(EGOA_re_LL$yrs,EGOA_re_LL$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,40000))
polygon(c(EGOA_re_LL$yrs,sort(EGOA_re_LL$yrs,decreasing=TRUE)),c(EGOA_re_LL$UCI,rev(EGOA_re_LL$LCI)),col="light grey",border=NA)
lines(EGOA_re_LL$yrs,EGOA_re_LL$biomA,lwd=3,col="orangered")
arrows(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est+1.96*EGOA_srv_LL$srv_est*EGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est-1.96*EGOA_srv_LL$srv_est*EGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,37770,"EGOA",cex=1.5)

axis(1,cex.axis=1.25)

#============================================================================================
#============= Plot Re fit to LL by region - same scale (for SAFE)
#============================================================================================

dev.new(width=9,height=7)
layout(matrix(c(0,1,2,3,0),5,1,byrow = TRUE),heights=c(0.1,1,1,1,0.2))
par(mar=c(0.1,7,0.1,1))

x<-plot(WGOA_re_LL$yrs,WGOA_re_LL$biomA,type="l",lwd=3,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,20000))
polygon(c(WGOA_re_LL$yrs,sort(WGOA_re_LL$yrs,decreasing=TRUE)),c(WGOA_re_LL$UCI,rev(WGOA_re_LL$LCI)),col="light grey",border=NA)
lines(WGOA_re_LL$yrs,WGOA_re_LL$biomA,lwd=3)
arrows(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est+1.96*WGOA_srv_LL$srv_est*WGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est-1.96*WGOA_srv_LL$srv_est*WGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(WGOA_srv_LL$yrs_srv,WGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,17770,"WGOA",cex=1.5)


x<-plot(CGOA_re_LL$yrs,CGOA_re_LL$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,40000))
polygon(c(CGOA_re_LL$yrs,sort(CGOA_re_LL$yrs,decreasing=TRUE)),c(CGOA_re_LL$UCI,rev(CGOA_re_LL$LCI)),col="light grey",border=NA)
lines(CGOA_re_LL$yrs,CGOA_re_LL$biomA,lwd=3)
arrows(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est+1.96*CGOA_srv_LL$srv_est*CGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est-1.96*CGOA_srv_LL$srv_est*CGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(CGOA_srv_LL$yrs_srv,CGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,37770,"CGOA",cex=1.5)
mtext("Relative Population Weight (RPW)",side=2,line=4.5,cex=1.25)

x<-plot(EGOA_re_LL$yrs,EGOA_re_LL$biomA,type="l",lwd=2,las=2,cex.axis=1.25,xaxt="n",ylab="",ylim=c(0,30000))
polygon(c(EGOA_re_LL$yrs,sort(EGOA_re_LL$yrs,decreasing=TRUE)),c(EGOA_re_LL$UCI,rev(EGOA_re_LL$LCI)),col="light grey",border=NA)
lines(EGOA_re_LL$yrs,EGOA_re_LL$biomA,lwd=3)
arrows(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est+1.96*EGOA_srv_LL$srv_est*EGOA_srv_LL$srv_sd,angle=90,length=0.05)
arrows(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est-1.96*EGOA_srv_LL$srv_est*EGOA_srv_LL$srv_sd,angle=90,length=0.05)
points(EGOA_srv_LL$yrs_srv,EGOA_srv_LL$srv_est,pch=16,cex=1.75,col="darkorange4",bg="white")
text(1985,27770,"EGOA",cex=1.5)

axis(1,cex.axis=1.25)

#============================================================================================
#============= Plot apportionment
#============================================================================================

Apport<-read.csv(paste(pathD,"/Apport.csv",sep=""))
Apport_17<-read.csv(paste(pathD,"/Apport_17.csv",sep=""))


dev.new(width=9,height=7)

x<-plot(Apport$yrs[1:35],Apport_17$WGOA,type="l",lwd=2,las=2,xaxt="n",ylab="",xlab="",ylim=c(0,0.5))
lines(Apport$yrs[1:35],Apport$WGOA[1:35],lwd=3,col="orangered")
text(Apport$yrs[35],Apport$WGOA[35],labels="W",pos=4)

lines(Apport$yrs[1:35],Apport_17$CGOA,lwd=2,lty=2)
lines(Apport$yrs[1:35],Apport$CGOA[1:35],lwd=3,col="orangered",lty=2)
text(Apport$yrs[35],Apport$CGOA[35],labels="C",pos=4)

lines(Apport$yrs[1:35],Apport_17$EGOA,lwd=2,lty=3)
lines(Apport$yrs[1:35],Apport$EGOA[1:35],lwd=3,col="orangered",lty=3)
text(Apport$yrs[35],Apport$EGOA[35],labels="E",pos=4)

axis(1)
mtext("Apportionment",side=2,line=2.5,cex=1.25)
legend(1984,0.1,legend=c("2015.1","2018.1"),col=c("black","orangered"),lty=c(1,1),lwd=c(2,3))
