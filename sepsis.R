library(dplyr)
library(tidyverse)
library(xlsx)
library(readxl)
library(openxlsx)


read_pre <- function(file){
  edcform = read_csv(file) 
  edcform %>%  rename_with(~ paste0(edcform$DataPageName[1],"_", .), -c("project","studyid", "Subject", "SDVTier","Site","SiteNumber","SiteGroup"))
      
}

####EDC
data_join <- list.files(path = "./EDC_pro/", # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_pre) %>%                              # Store all files in list
  reduce(full_join, by = "Subject")                      # Full-join data sets into one data set 

###instrument
ins <- read_csv("Interim1_instruments.csv")

mdw_filter = ins %>% 
select(RDFilename,matches("MDW")) %>% 
filter(Diff_MDW_Flag_NonNumericFlag=="NULL"& Diff_MDW_SingleCharParameterFlag=="NULL") %>%
separate(RDFilename, into = c('Time','标本编号'), sep = '_') 
# left_join(Patient_2023_07_04_18h09m_IVD, by = "标本编号") %>%


###clinical label

label <- read_excel("Interim_label1.xlsx", skip = 1)

label1 = label[!is.na(label[,2]),1:2]
colnames(label1)=c("Subject","Label")

label2 = label[!is.na(label[,4]),3:4]
colnames(label2)=c("Subject","Label")

label3 = label[!is.na(label[,6]),5:6]
colnames(label3)=c("Subject","Label")

Sepsis_Label <- rbind( rbind(label1 ,  label2), label3) 


res_t = data_join %>% rename("标本编号"="Patient Information_DXHSMP") %>%
       # unite(标本编号, c("StudyEnvSiteNumber", "ENROLL_NUM"),sep = "") %>% 
       right_join(mdw_filter, by = "标本编号") %>% 
       left_join(Sepsis_Label, by = "Subject") 




###generate final table
# res_t$Duplicate <- duplicated(res_t$Subject)
res = res_t %>% drop_na(Subject) %>% 
      rename(CBCADAT="Blood Collection for the Study_CBCADAT") %>%
      mutate(CBCADAT=format(strptime(CBCADAT, "%d/%b/%Y %H:%M"),format='%Y-%m-%d %H:%M')) %>% 
      mutate(Time=format(as_datetime(Time),format='%Y-%m-%d %H:%M')) %>% 
      mutate(Time_Check = ifelse(CBCADAT==Time, "Yes", "No")) %>%
      group_by(Subject) %>% mutate(freq=n()) %>%
      filter(!(freq!=1 & CBCADAT!=Time)) %>%
      select(Subject,标本编号,Time_Check,CBCADAT,Time,Enrollment_ENROLLYN,Label,matches("Site")[1],Diff_MDW_Value, "CEC Adjudicator 1_SFDIAGA","CEC Adjudicator 1_FSDIAGARB","CEC Adjudicator 2_SFDIAGA","CEC Adjudicator 2_FSDIAGARB","CEC Arbitrator_SFDIAGA","CEC Arbitrator_FSDIAGARB") 

colnames(res)=c("Subject", "标本编号","仪器分析时间检查","全血细胞分类计数分析日期时间","仪器真实分析时间","入组","剔除标签","Site","MDW", "Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3")

write.xlsx(res,  "Sepsis_STAT.xlsx",  colNames = TRUE)
