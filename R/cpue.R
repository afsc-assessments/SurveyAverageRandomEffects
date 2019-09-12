#============================================================================================
#============= Get Working Directories
#============================================================================================
source("read-admb.R")
library(ggplot2)
library(ggthemes)
library(tidyr)
library(dplyr)
library(scales)


titles <-c( "CPUE to survey weight 0:1 ", "CPUE to survey weight 1:2 ", "CPUE to survey weight 1:1 ", "CPUE to survey weight 2:1 ", "CPUE to survey weight 1:0 ")
area <- c("Eastern","Central","Western")
fn <- c( "cpue_srv_0_1.rep", "cpue_srv_1_2.rep", "cpue_srv_1_1.rep", "cpue_srv_2_1.rep", "cpue_srv_1_0.rep")
df1<-NULL
j=2
for (j in 1:5){
	ml <- read_rep(fn[j])
	for (i in 1:3){
		df2 <- data.frame(Year=ml$yrs, Area=area[i],
		source = "RE model", biom = ml$biomA[,i], uci  = ml$UCI[,i], lci  = ml$LCI[,i], case = j) 
		df1 <- rbind(df1,df2)
	}
	df1 <- rbind(df1,data.frame(Year=ml$yrs,Area="Combined",
	source = "RE model", biom = ml$biom_TOT, uci  = ml$biom_TOT_UCI, lci  = ml$biom_TOT_LCI, case=j) )

  p1 <- df1 %>% filter(case==j, Area != "Combined") %>% ggplot(aes(x=Year, y=biom,ymax=uci,ymin=lci ,fill=Area)) + geom_line() +
     geom_ribbon(alpha=.7) + theme_few( base_size=16) + facet_grid(Area~.,scales="free") + ggtitle(titles[j])
     ylab("Biomass")
  ggsave(paste0("mod_",j,".png"))
  p2 <- df1 %>% filter(case==j, Area != "Combined") %>% group_by(Year) %>% mutate(proportion=biom/sum(biom)) %>% 
     ggplot(aes(x=Year, y=proportion,fill=Area)) + geom_area() + ggtitle(titles[j]) +
     theme_few( base_size=16) + ylab("Biomass")
  ggsave(paste0("app_",j,".png"))
}

# Make a table from big dataframe, rows are models, columns are area-specific % for 2019
tt <- df1 %>% filter(Year==2019, Area != "Combined") %>% group_by(Year,case) %>% transmute(Area,proportion=biom/sum(biom)) 
dd <- spread(tt,Area,proportion) 
library(flextable)
kable
ggplot(aes(x=Year, y=biom,ymax=uci,ymin=lci ,fill=Area)) + geom_line() +

#"yrs_srv"         "srv_est_TOT"     "yrs_srv_LL"      "srv_est_TOT_LL"  "log_q_LL"        "yrs"             "biom_TOT"       
#"SD_biom_TOT"     "biom_TOT_UCI"    "biom_TOT_LCI"    "biom_TOT_LL"     "SD_biom_TOT_LL"  "biom_TOT_UCI_LL" "srv_est"        
#"srv_est_LL"      "LCI"             "biomA"           "UCI"             "biomsd"          "biomsd.sd"      
