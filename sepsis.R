library(dplyr)
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(skimr)

read_pre <- function(file){
  print(file)
  edcform = read_csv(file) 
  edcform %>% select(-c("projectid","project","studyid","environmentName","subjectId","StudySiteId","SDVTier","siteid","Site","SiteNumber","SiteGroup","instanceId","InstanceName","InstanceRepeatNumber","folderid","Folder","FolderName","FolderSeq","TargetDays","DataPageId","PageRepeatNumber","RecordDate","RecordId","RecordPosition","MinCreated","MaxUpdated","SaveTS","StudyEnvSiteNumber")) %>%
  rename_with(~ paste0(edcform$DataPageName[1],"_", .), -c("Subject","DataPageName"))
}

####EDC
data_join <- list.files(path = "./EDC/", # Identify all CSV files
                       pattern = "*.CSV", full.names = TRUE) %>% 
  lapply(read_pre) %>%                              # Store all files in list
  reduce(full_join, by = "Subject")                      # Full-join data sets into one data set 

###instrument
ins <- read_csv("Interim2_instruments.csv")

mdw_filter = ins %>% 
select(RDFilename,matches("MDW"),matches("Ly_DC"),matches("Ly_Op_Mean")) %>% 
filter(Diff_MDW_Flag_NonNumericFlag=="NULL"& Diff_MDW_SingleCharParameterFlag=="NULL") %>%
separate(RDFilename, into = c('Time','标本编号'), sep = '_') 
# left_join(Patient_2023_07_04_18h09m_IVD, by = "标本编号") %>%

ins %>% select(matches("Ly"))  %>% mutate_if(is.character, as.numeric) %>% skimr::skim()



###clinical label

label <- read_excel("Interim2_instruments.xlsx", skip = 1)

label1 = label[!is.na(label[,2]),1:3]
colnames(label1)=c("Subject","Label","Batch")

label2 = label[!is.na(label[,4]),4:6]
colnames(label2)=c("Subject","Label","Batch")

label3 = label[!is.na(label[,6]),7:9]
colnames(label3)=c("Subject","Label","Batch")

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
      #compute Lymph_index
      mutate_at(c('Ly_DC_Mean', 'Ly_DC_SD', 'Ly_Op_Mean'), as.numeric) %>% mutate(Lymph_index =Ly_DC_Mean*Ly_DC_SD/Ly_Op_Mean) %>% 
      mutate(Label2 = ifelse(Lymph_index > 11.68, "可能感染", NA)) %>% 
      select(Subject,标本编号,Time_Check,CBCADAT,Time,Enrollment_ENROLLYN,Label,Label2,Lymph_index,matches("Site")[1],Diff_MDW_Value,"CEC Adjudicator 1_SFDIAGA","CEC Adjudicator 1_FSDIAGARB","CEC Adjudicator 2_SFDIAGA","CEC Adjudicator 2_FSDIAGARB","CEC Arbitrator_SFDIAGA","CEC Arbitrator_FSDIAGARB") 

colnames(res)=c("Subject", "标本编号","仪器分析时间检查","全血细胞分类计数分析日期时间","仪器真实分析时间","入组","剔除标签","新冠感染","Lymph_index","Site","MDW", "Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3")

write.xlsx(res,  "Sepsis_STAT_Covid.xlsx",  colNames = TRUE)


######explore

Sepsis_type=c("Septic shock (severe sepsis + hypotension)","Severe sepsis (organ dysfunction or tissue hypoperfusion)","Sepsis (2 SIRS + infection)")

res2 = res %>% filter(is.na(剔除标签)) %>%
       mutate_at(c('MDW'), as.numeric) %>% 
       drop_na(Adjudicator1_Sepsis2) %>% drop_na(Adjudicator2_Sepsis2) %>%
       mutate(Sepsis2=ifelse(Adjudicator1_Sepsis2 == Adjudicator2_Sepsis2, Adjudicator1_Sepsis2, Arbitrator_Sepsis2))%>% 
       drop_na(Sepsis2) %>% 
       mutate(Sepsis2_Sum = ifelse(Sepsis2 %in% Sepsis_type, "Sepsis" , Sepsis2))

my_comparisons <- list( c("Infection without sepsis", "Sepsis"), c("Infection without sepsis", "SIRS (>= 2 SIRS criteria)"), c("Non-SIRS/non-infection (control case)","Infection without sepsis") )


outpdf=paste("Sepsis","_cor.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("Non-SIRS/non-infection (control case)","SIRS (>= 2 SIRS criteria)", "Infection without sepsis", "Sepsis")),y=Lymph_index)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$Lymph_index), linetype = 2) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")


ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("Non-SIRS/non-infection (control case)","SIRS (>= 2 SIRS criteria)", "Infection without sepsis", "Sepsis")),y=MDW)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$MDW), linetype = 2) +
stat_compare_means( ref.group = ".all.",method = "wilcox.test", label = "p.signif")


ggplot(res2, aes(x = MDW, y = Lymph_index)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "blue", alpha = 0.3) +
  theme_bw()+stat_cor(data=res2, method = "spearman")

ggplot(res2, aes(x = MDW, y = Lymph_index, color=Sepsis2_Sum)) + geom_point() + stat_smooth(method = 'lm')

dev.off()


res3 =  res2 %>% filter(Sepsis2_Sum== "Non-SIRS/non-infection (control case)") %>%
        mutate(MDW_Diag = ifelse(MDW > 20.5, "假阳", NA)) %>% 
        mutate(Ly_Diag = ifelse(Lymph_index > 11.68, "病毒感染", NA)) %>%
        select(MDW,MDW_Diag,Lymph_index,Ly_Diag,Sepsis2_Sum,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 
        
write.xlsx(res3,  "Non-SIRS_screen.xlsx",  colNames = TRUE)


res4 =  res2 %>% filter(Sepsis2_Sum== "SIRS (>= 2 SIRS criteria)") %>%
        mutate(MDW_Diag = ifelse(MDW > 20.5, "假阳", NA)) %>% 
        mutate(Ly_Diag = ifelse(Lymph_index > 11.68, "病毒感染", NA)) %>%
        select(MDW,MDW_Diag,Lymph_index,Ly_Diag,Sepsis2_Sum,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 
        
write.xlsx(res4,  "SIRS_screen.xlsx",  colNames = TRUE)