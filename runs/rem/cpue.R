#============================================================================================
#============= Get Working Directories
#============================================================================================
getwd()
source("../../R/read-admb.R")
source("../../R/read-admb2.R")
library(ggplot2)
library(ggthemes)
library(tidyr)
library(dplyr)
library(scales)

write_dat <- function(data,fn) {
  sink(fn)
  for (i in 1:length(data)){
    cat(paste0("#",names(data[i]),"\n"))
  	d <- data[[i]]
    if(!is.null(dim(d) )){  
    	write.table(d,row.names=F,col.names=F) 
    } else{
    	cat(paste0(d,"\n"))
    }
  }
  sink()
}
d1 <-	read_dat("atka.dat"); 
wt  = c(1,0.666,0.5,0.333,0.001)
wt2 = c(1,0.666,0.5,0.333,0.001)
d1
runmod <- function(i,j){
	d1$srv_wt  <- wt[i]
  d1$srv2_wt <-wt2[j]
  write_dat(d1,"rem.dat")
  system("./rem -nox")
  file.copy("rwout.rep",paste0("wt_",i,"_",j,".rep"),overwrite=TRUE)
}
i=3
for (i in 1:5){
  j=(6-i)
  runmod(i,j)
} 

titles <-c( "CPUE weight 0 ", "CPUE weight 1:2 ", "CPUE weight 1:1 ", "CPUE weight 2:1 ", "CPUE weight 1:0 ")
area <- c("Eastern","Central","Western")
fn <- c( "wt_1_5.rep", "wt_2_4.rep", "wt_3_3.rep", "wt_4_2.rep", "wt_5_1.rep")
df1<-NULL
j=4
for (j in 1:5){
	ml <- read_rep(fn[j])
	for (i in 1:3){
		df2 <- data.frame(Year=ml$yrs, Area=area[i],
		source = "RE model", biom = ml$biomA[,i], uci  = ml$UCI[,i], lci  = ml$LCI[,i], case = j) 
		df3 <- data.frame(Year=ml$yrs_srv, Area=area[i], source = "Bottom trawl survey", biom = ml$srv_est[,i], uci  = NA, lci  = NA, case = j) 
		df4 <- data.frame(Year=ml$yrs_srv2, Area=area[i], source = "Fishery CPUE", biom = ml$srv_est_srv2[,i]/exp(ml$log_q_srv2[i]), uci  = NA, lci  = NA, case = j) 
		df1 <- rbind(df1,df2,df3,df4)
	}
	df1 <- rbind(df1,data.frame(Year=ml$yrs,Area="Combined",
	source = "RE model", biom = ml$biom_TOT, uci  = ml$biom_TOT_UCI, lci  = ml$biom_TOT_LCI, case=j) )

  p1 <- df1 %>% filter(source=="RE model",case==j, Area != "Combined") %>% 
     ggplot(aes(x=Year, y=biom,ymax=uci,ymin=lci )) + geom_line() +
     geom_ribbon(alpha=.3,fill='gold') + theme_few( base_size=16) + facet_grid(Area~.,scales="free") + ggtitle(titles[j]) +
     ylab("Biomass") + expand_limits(y = 0)
  dd <- df1 %>% filter(source!="RE model",case==j, Area != "Combined") 
#  d2 <- df1 %>% filter(source=="Fishery CPUE",case==j, Area != "Combined") 
  p1 <- p1 + geom_point(data=dd,aes(color=source),size=4)
  p1
  ggsave(paste0("mod_",j,".png"))
  head(df1)

  tmp <- df1 %>% filter(source=="RE model",case==j, Area != "Combined") %>% group_by(Year) %>% mutate(proportion=biom/sum(biom)) 
  head(tmp)
  p2 <- ggplot(tmp,aes(x=Year, y=proportion,fill=Area)) + geom_area() + ggtitle(titles[j]) +
     theme_few( base_size=16) + ylab("Biomass")#+ expand_limits(y = 0) 
    p2
  ggsave(paste0("app_",j,".png"))
}
head(df1)
unique(df1$source)

# Plot fits to data
		df3 <- data.frame(Year=ml$yrs_srv, Area=area[i],
		source = "Survey", biom = ml$srv_est[,i], uci  = ml$UCI[,i], lci  = ml$LCI[,i], case = j) 
		df1 <- rbind(df1,df2)

head(tt)
# Make a table from big dataframe, rows are models, columns are area-specific % for 2019
tt <- df1 %>% filter(source=="RE model",Year==2019, Area != "Combined") %>% group_by(Year,case) %>% transmute(Area,proportion=biom/sum(biom)) 
dd <- spread(tt,Area,proportion) 
dd
library(flextable)
kable
ggplot(aes(x=Year, y=biom,ymax=uci,ymin=lci ,fill=Area)) + geom_line() +

#"yrs_srv"         "srv_est_TOT"     "yrs_srv_LL"      "srv_est_TOT_LL"  "log_q_LL"        "yrs"             "biom_TOT"       
#"SD_biom_TOT"     "biom_TOT_UCI"    "biom_TOT_LCI"    "biom_TOT_LL"     "SD_biom_TOT_LL"  "biom_TOT_UCI_LL" "srv_est"        
#"srv_est_LL"      "LCI"             "biomA"           "UCI"             "biomsd"          "biomsd.sd"      
