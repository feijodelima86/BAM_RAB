library(readr)
library(tidyr)
library(zoo)
library(dplyr)

FINAL_PARTIAL_REDUX <- read.csv("0_data/icp_ms_partial/FINAL_PARTIAL_REDUX.csv")

names(FINAL_PARTIAL_REDUX)

df2<-FINAL_PARTIAL_REDUX[,(c(1:3,6:15))]


df3 <-aggregate(x = df2[,colnames(df2) != c("Date","Site","Fraction")],          
                  by = list(df2$Date,df2$Site,df2$Fraction),
                  FUN = mean,
                  na.rm = F)

write.csv(df3, "2_incremental/water1.csv")


SIZE_FRACTIONS <- read.csv("2_incremental/water2.csv")

SIZE_FRACTIONS[SIZE_FRACTIONS < 0] <- 0  

names(SIZE_FRACTIONS)

n=1

PS<-SIZE_FRACTIONS[,n]  -(SIZE_FRACTIONS[,n+13])
PC<-SIZE_FRACTIONS[,n+13] -(SIZE_FRACTIONS[,n+26])
PD<-SIZE_FRACTIONS[,n+26]

S1<-data.frame(SIZE_FRACTIONS[,c(2:3)])

for(i in 5:14){
  S1[,i-2]<- SIZE_FRACTIONS[,i] - SIZE_FRACTIONS[,i+13]
  S1[,i+8]<- SIZE_FRACTIONS[,i+13] - SIZE_FRACTIONS[,i+26]
  S1[,i+18]<- SIZE_FRACTIONS[,i+26]
  }

S1[,i]<- SIZE_FRACTIONS[,i+13] - SIZE_FRACTIONS[,i+26]

var_names<-c("P","Ca","Fe","Cu","Zn","As","Se","Mo","Cd","Pb")

colnames(S1)[c(3:12)]<-paste("SS", var_names, sep=" ")
colnames(S1)[c(13:22)]<-paste("CL", var_names, sep=" ")
colnames(S1)[c(23:32)]<-paste("TD", var_names, sep=" ")

S1[S1 < 0] <- 0 

write.csv(S1, "2_incremental/water3.csv")




