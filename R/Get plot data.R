#============================================================================================
#============= Get Working Directories
#============================================================================================

path<-getwd()

pathM<-paste(path,"",sep="")
pathD<-paste(path,"/Plot data",sep="")

path
#============================================================================================
#============= Read in results data file, set up results matrices
#============================================================================================

RWOUT<-readLines(paste(pathM,"/RWOUT.rep",sep=""),warn=FALSE)
(RWOUT)

yrs<-as.numeric(strsplit(RWOUT[grep("yrs",RWOUT)[3]+1]," ")[[1]][2:length(strsplit(RWOUT[grep("yrs",RWOUT)[3]+1]," ")[[1]])])
yrs_srv<-as.numeric(strsplit(RWOUT[grep("yrs_srv",RWOUT)[1]+1]," ")[[1]][2:length(strsplit(RWOUT[grep("yrs_srv",RWOUT)[1]+1]," ")[[1]])])
yrs_cpue<-as.numeric(strsplit(RWOUT[grep("yrs_srv_LL",RWOUT)[1]+1]," ")[[1]][2:length(strsplit(RWOUT[grep("yrs_srv_LL",RWOUT)[1]+1]," ")[[1]])])
q_cpue<-exp(as.numeric(strsplit(RWOUT[grep("log_q_LL",RWOUT)+1]," ")[[1]][2:length(strsplit(RWOUT[grep("log_q_LL",RWOUT)+1]," ")[[1]])]))
q_cpue
yrs_cpue

# RE estimates
df1 <- data.frame(yrs=yrs,Area="Combined",
	source = "RE model",
	biom = as.numeric(strsplit(RWOUT[grep("biom_TOT",RWOUT)[1]+1]," ")[[1]][2:(length(yrs)+1)]),
  uci  = as.numeric(strsplit(RWOUT[grep("biom_TOT_UCI",RWOUT)[1]+1]," ")[[1]][2:(length(yrs)+1)]),
  lci  = as.numeric(strsplit(RWOUT[grep("biom_TOT_LCI",RWOUT)[1]+1]," ")[[1]][2:(length(yrs)+1)]) )

df2 <- data.frame(yrs=yrs,Area="Combined",
	source = "CPUE",
  biom   = as.numeric(strsplit(RWOUT[grep("biom_TOT_LL",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]),
  uci    = as.numeric(strsplit(RWOUT[grep("biom_TOT_UCI_LL",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]),
  lci    = as.numeric(strsplit(RWOUT[grep("biom_TOT_LCI_LL",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]) )

df3 <- data.frame(yrs=yrs,Area="Combined",
	source = "Survey",
  biom   = as.numeric(strsplit(RWOUT[grep("biom_TOT_srv",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]),
  uci    = as.numeric(strsplit(RWOUT[grep("biom_TOT_UCI_srv",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]),
  lci    = as.numeric(strsplit(RWOUT[grep("biom_TOT_LCI_srv",RWOUT)+1]," ")[[1]][2:(length(yrs)+1)]) )

EAI_re<-matrix(nrow=length(yrs),ncol=4)
colnames(EAI_re)<-c("yrs","biomA","UCI","LCI")
EAI_re[,1]<-yrs
EAI_re_cpue<-matrix(nrow=length(yrs),ncol=4)
colnames(EAI_re_cpue)<-c("yrs","biomA","UCI","LCI")
EAI_re_cpue[,1]<-yrs
EAI_srv<-matrix(nrow=length(yrs_srv),ncol=3)
colnames(EAI_srv)<-c("yrs_srv","srv_est","srv_sd")
EAI_srv[,1]<-yrs_srv
EAI_cpue<-matrix(nrow=length(yrs_cpue),ncol=3)
colnames(EAI_cpue)<-c("yrs_srv","srv_est","srv_sd")
EAI_cpue[,1]<-yrs_cpue

Apport<-matrix(nrow=length(yrs),ncol=4)
colnames(Apport)<-c("yrs","WAI","CAI","EAI")
Apport[,1]<-yrs
ls()
#============================================================================================
#============= Pull AI-wide results
#============================================================================================

# Trawl survey data
AI_srv[,2]<-as.numeric(strsplit(RWOUT[grep("srv_est_TOT",RWOUT)+1]," ")[[1]][2:length(strsplit(RWOUT[grep("srv_est_TOT",RWOUT)+1]," ")[[1]])])

reg_strat<-length(as.numeric(strsplit(RWOUT[grep("srv_est",RWOUT)[3]+1]," ")[[1]]))-1
num_strat<-reg_strat/3
srv_est<-matrix(nrow=length(yrs_srv),ncol=reg_strat)
srv_var<-matrix(nrow=length(yrs_srv),ncol=reg_strat)

for(y in 1:length(yrs_srv)){
srv_est[y,]<-as.numeric(strsplit(RWOUT[(grep("srv_est",RWOUT)[3]+1):(grep("srv_est",RWOUT)[3]+length(yrs_srv))]," ")[[y]][2:(reg_strat+1)])
srv_est[y,which(srv_est[y,]<0)]<-0
srv_sd<-as.numeric(strsplit(RWOUT[(grep("srv_sd",RWOUT)[1]+1):(grep("srv_sd",RWOUT)[1]+length(yrs_srv))]," ")[[y]][2:(reg_strat+1)])
srv_sd_sd<-(srv_est[y,]+0.0001)*sqrt(exp(srv_sd^2)-1)
srv_sd_sd[which(srv_sd_sd<=0.1)]<-0
srv_var[y,]<-srv_sd_sd^2}

srv_var_TOT<-rowSums(srv_var)
AI_srv[,3]<-sqrt(srv_var_TOT)/AI_srv[,2]

# CPUE survey data
AI_cpue[,2]<-as.numeric(strsplit(RWOUT[grep("srv_est_TOT_cpue",RWOUT)+1]," ")[[1]][2:length(strsplit(RWOUT[grep("srv_est_TOT_cpue",RWOUT)+1]," ")[[1]])])

reg_strat_cpue<-length(as.numeric(strsplit(RWOUT[grep("srv_est_cpue",RWOUT)+1]," ")[[1]]))-1
num_strat_cpue<-reg_strat_cpue/3
srv_est_cpue<-matrix(nrow=length(yrs_cpue),ncol=reg_strat_cpue)
srv_var_cpue<-matrix(nrow=length(yrs_cpue),ncol=reg_strat_cpue)

for(y in 1:length(yrs_cpue)){
srv_est_cpue[y,]<-as.numeric(strsplit(RWOUT[(grep("srv_est_cpue",RWOUT)+1):(grep("srv_est_cpue",RWOUT)+length(yrs_cpue))]," ")[[y]][2:(reg_strat_cpue+1)])
srv_est_cpue[y,which(srv_est_cpue[y,]<0)]<-0
srv_sd_cpue<-as.numeric(strsplit(RWOUT[(grep("srv_sd_cpue",RWOUT)[1]+1):(grep("srv_sd_cpue",RWOUT)[1]+length(yrs_cpue))]," ")[[y]][2:(reg_strat_cpue+1)])
srv_sd_sd_cpue<-(srv_est_cpue[y,]+0.0001)*sqrt(exp(srv_sd_cpue^2)-1)
srv_sd_sd_cpue[which(srv_sd_sd_cpue<=0.1)]<-0
srv_var_cpue[y,]<-srv_sd_sd_cpue^2}

srv_var_TOT_cpue<-rowSums(srv_var_cpue)
AI_cpue[,3]<-sqrt(srv_var_TOT_cpue)/AI_cpue[,2]


#============================================================================================
#============= Pull W/C/EAI results
#============================================================================================

# Trawl survey data
CAI_srv[,2]<-srv_est[,1]
CAI_srv[,3]<-sqrt(srv_var[,1])/CAI_srv[,2]
EAI_srv[,2]<-srv_est[,2]
EAI_srv[,3]<-sqrt(srv_var[,2])/EAI_srv[,2]
EAI_srv<-EAI_srv[-which(EAI_srv[,2]==0),]
WAI_srv[,2]<-srv_est[,3]
WAI_srv[,3]<-sqrt(srv_var[,3])/WAI_srv[,2]

# cpue survey data
CAI_cpue[,2]<-srv_est_cpue[,1]
CAI_cpue[,3]<-sqrt(srv_var_cpue[,1])/CAI_cpue[,2]
EAI_cpue[,2]<-srv_est_cpue[,2]
EAI_cpue[,3]<-sqrt(srv_var_cpue[,2])/EAI_cpue[,2]
WAI_cpue[,2]<-srv_est_cpue[,3]
WAI_cpue[,3]<-sqrt(srv_var_cpue[,3])/WAI_cpue[,2]

# RE estimates
biomA<-matrix(nrow=length(yrs),ncol=reg_strat)
biomsd<-matrix(nrow=length(yrs),ncol=reg_strat)
biomsd_sd<-matrix(nrow=length(yrs),ncol=reg_strat)

for(y in 1:length(yrs)){
biomA[y,]<-as.numeric(strsplit(RWOUT[(grep("biomA",RWOUT)+1):(grep("biomA",RWOUT)+length(yrs))]," ")[[y]][2:(reg_strat+1)])
biomsd[y,]<-as.numeric(strsplit(RWOUT[(grep("biomsd",RWOUT)[1]+1):(grep("biomsd",RWOUT)[1]+length(yrs))]," ")[[y]][2:(reg_strat+1)])
biomsd_sd[y,]<-as.numeric(strsplit(RWOUT[(grep("biomsd.sd",RWOUT)+1):(grep("biomsd.sd",RWOUT)+length(yrs))]," ")[[y]][2:(reg_strat+1)])}

# CAI
numer<-exp(2*biomsd[,1]+biomsd_sd[,1]^2)*(exp(biomsd_sd[,1]^2)-1)
denom<-exp(biomsd[,1]+0.5*biomsd_sd[,1]^2)^2
SD_biom_TOT<-sqrt(log(numer/denom+1))
CAI_re[,2]<-biomA[,1]
CAI_re[,3]<-exp(log(CAI_re[,2])+1.96*SD_biom_TOT)
CAI_re[,4]<-exp(log(CAI_re[,2])-1.96*SD_biom_TOT)
numer_cpue<-exp(2*log(q_cpue[1]*exp(biomsd[,1]))+biomsd_sd[,1]^2)*(exp(biomsd_sd[,1]^2)-1)
denom_cpue<-exp(log(q_cpue[1]*exp(biomsd[,1]))+0.5*biomsd_sd[,1]^2)^2
SD_biom_TOT_cpue<-sqrt(log(numer_cpue/denom_cpue+1))
CAI_re_cpue[,2]<-q_cpue[1]*biomA[,1]
CAI_re_cpue[,3]<-exp(log(q_cpue[1]*CAI_re[,2])+1.96*SD_biom_TOT_cpue)
CAI_re_cpue[,4]<-exp(log(q_cpue[1]*CAI_re[,2])-1.96*SD_biom_TOT_cpue)


# EAI
numer<-exp(2*biomsd[,2]+biomsd_sd[,2]^2)*(exp(biomsd_sd[,2]^2)-1)
denom<-exp(biomsd[,2]+0.5*biomsd_sd[,2]^2)^2
SD_biom_TOT<-sqrt(log(numer/denom+1))
EAI_re[,2]<-biomA[,2]
EAI_re[,3]<-exp(log(EAI_re[,2])+1.96*SD_biom_TOT)
EAI_re[,4]<-exp(log(EAI_re[,2])-1.96*SD_biom_TOT)
numer_cpue<-exp(2*log(q_cpue[2]*exp(biomsd[,2]))+biomsd_sd[,2]^2)*(exp(biomsd_sd[,2]^2)-1)
denom_cpue<-exp(log(q_cpue[2]*exp(biomsd[,2]))+0.5*biomsd_sd[,2]^2)^2
SD_biom_TOT_cpue<-sqrt(log(numer_cpue/denom_cpue+1))
EAI_re_cpue[,2]<-q_cpue[2]*biomA[,2]
EAI_re_cpue[,3]<-exp(log(q_cpue[2]*EAI_re[,2])+1.96*SD_biom_TOT_cpue)
EAI_re_cpue[,4]<-exp(log(q_cpue[2]*EAI_re[,2])-1.96*SD_biom_TOT_cpue)

# WAI
numer<-exp(2*biomsd[,3]+biomsd_sd[,3]^2)*(exp(biomsd_sd[,3]^2)-1)
denom<-exp(biomsd[,3]+0.5*biomsd_sd[,3]^2)^2
SD_biom_TOT<-sqrt(log(numer/denom+1))
WAI_re[,2]<-biomA[,3]
WAI_re[,3]<-exp(log(WAI_re[,2])+1.96*SD_biom_TOT)
WAI_re[,4]<-exp(log(WAI_re[,2])-1.96*SD_biom_TOT)
numer_cpue<-exp(2*log(q_cpue[3]*exp(biomsd[,3]))+biomsd_sd[,3]^2)*(exp(biomsd_sd[,3]^2)-1)
denom_cpue<-exp(log(q_cpue[3]*exp(biomsd[,3]))+0.5*biomsd_sd[,3]^2)^2
SD_biom_TOT_cpue<-sqrt(log(numer_cpue/denom_cpue+1))
WAI_re_cpue[,2]<-q_cpue[3]*biomA[,3]
WAI_re_cpue[,3]<-exp(log(q_cpue[3]*WAI_re[,2])+1.96*SD_biom_TOT_cpue)
WAI_re_cpue[,4]<-exp(log(q_cpue[3]*WAI_re[,2])-1.96*SD_biom_TOT_cpue)

# Compute apportionment
reg_biom<-cbind(WAI_re[,2],CAI_re[,2],EAI_re[,2])
for(y in 1:length(yrs)){
Apport[y,2:4]<-reg_biom[y,]/sum(reg_biom[y,])}


#============================================================================================
#============= Write plot data
#============================================================================================

write.csv(CAI_re,paste(pathD,"/CAI_re.csv",sep=""))
write.csv(CAI_srv,paste(pathD,"/CAI_srv.csv",sep=""))
write.csv(CAI_re_cpue,paste(pathD,"/CAI_re_cpue.csv",sep=""))
write.csv(CAI_cpue,paste(pathD,"/CAI_cpue",sep=""))
write.csv(WAI_re,paste(pathD,"/WAI_re.csv",sep=""))
write.csv(WAI_srv,paste(pathD,"/WAI_srv.csv",sep=""))
write.csv(WAI_re_cpue,paste(pathD,"/WAI_re_cpue.csv",sep=""))
write.csv(WAI_cpue,paste(pathD,"/WAI_cpue",sep=""))
write.csv(EAI_re,paste(pathD,"/EAI_re.csv",sep=""))
write.csv(EAI_srv,paste(pathD,"/EAI_srv.csv",sep=""))
write.csv(EAI_re_cpue,paste(pathD,"/EAI_re_cpue.csv",sep=""))
write.csv(EAI_cpue,paste(pathD,"/EAI_cpue",sep=""))
write.csv(AI_re,paste(pathD,"/AI_re.csv",sep=""))
write.csv(AI_srv,paste(pathD,"/AI_srv.csv",sep=""))
write.csv(AI_re_cpue,paste(pathD,"/AI_re_cpue.csv",sep=""))
write.csv(AI_cpue,paste(pathD,"/AI_cpue",sep=""))
write.csv(Apport,paste(pathD,"/Apport.csv",sep=""))

