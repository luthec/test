library(dplyr)
library(tidyverse)
library(xlsx)
library(readxl)
library(openxlsx)

####EDC
data_join <- list.files(path = "./EDC_pro/", # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Subject")                      # Full-join data sets into one data set 

###instrument
all_test <- read_csv("all_test.csv")
Patient_2023_07_04_18h09m_IVD <- read_csv("Patient_2023-07-04_18h09m_IVD.csv")


mdw_filter = all_test %>% 
select(RDFilename,matches("MDW")) %>% 
filter(Diff_MDW_Flag_NonNumericFlag=="NULL"& Diff_MDW_SingleCharParameterFlag=="NULL") %>%
separate(RDFilename, into = c('Time','标本编号'), sep = '_') %>%
left_join(Patient_2023_07_04_18h09m_IVD, by = "标本编号") %>%
select(Time,标本编号,matches("MDW"))


mdw_filter2 = all_test %>% 
select(RDFilename,matches("MDW")) %>% 
filter(Diff_MDW_Flag_NonNumericFlag=="NULL"& Diff_MDW_SingleCharParameterFlag=="NULL") %>%
separate(RDFilename, into = c('Time','标本编号'), sep = '_') %>%
select(Time,标本编号,"Diff_MDW_Value")



write.xlsx(mdw_filter2 ,  "Instrument_test2.xlsx",  col.names = TRUE)


res_t = data_join %>% select( Subject,DXHSMP, CBCADAT) %>% rename("标本编号"="DXHSMP") %>%
       # unite(标本编号, c("StudyEnvSiteNumber", "ENROLL_NUM"),sep = "") %>% 
       right_join(Patient_2023_07_04_18h09m_CPD, by = "标本编号") %>% 
       select(Subject,标本编号,CBCADAT,分析日期,分析时间)

write.xlsx(res_t,  "Sepsis_CBCADAT.xlsx",  col.names = TRUE)


Sepsis_MDW_raw <- read_excel("SepsisMDW.xlsx", skip = 1)

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
