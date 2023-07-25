library(dplyr)
library(tidyverse)
library(xlsx)
library(readxl)
library(openxlsx)

CHN_097_9730 <- read_csv("CHN-097_9730.csv")
CHN_097_9732 <- read_csv("CHN-097_9732.csv")
CHN_097_9728 <- read_csv("CHN-097_9728.csv")


Sepsis_MDW_raw <- read_excel("Sepsis_MDW.xlsx", skip = 1)

Sepsis_MDW1 = Sepsis_MDW_raw[!is.na(Sepsis_MDW_raw[,2]),1:2]
colnames(Sepsis_MDW1)=c("Subject","MDW")

Sepsis_MDW2 = Sepsis_MDW_raw[!is.na(Sepsis_MDW_raw[,6]),5:6]
colnames(Sepsis_MDW2)=c("Subject","MDW")

Sepsis_MDW3 = Sepsis_MDW_raw[!is.na(Sepsis_MDW_raw[,10]),9:10]
colnames(Sepsis_MDW3)=c("Subject","MDW")

Sepsis_MDW  <- rbind( rbind(Sepsis_MDW1,  Sepsis_MDW2), Sepsis_MDW3) 


res1 = full_join(Sepsis_MDW,CHN_097_9730,by="Subject")   

res2 = full_join(res1 ,CHN_097_9732,by="Subject")

res3 = full_join(res2 ,CHN_097_9728,by="Subject")

res = res3[!is.na(res3$project.x),] %>% select( Subject, MDW, Site.x, DataPageName.x, SFDIAGA.x, FSDIAGARB.x,DataPageName.y, SFDIAGA.y, FSDIAGARB.y,DataPageName,SFDIAGA, FSDIAGARB)

colnames(res)=c("Subject", "MDW", "Site","Adjudicator1", "Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2", "Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator", "Arbitrator_Sepsis2", "Arbitrator_Sepsis3")

write.xlsx(res,  "Sepsis_STAT.xlsx",  col.names = TRUE)
