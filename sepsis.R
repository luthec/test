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
  edcform = read_csv(file,show_col_types = FALSE) 
  print(dim(edcform))
  edcform %>% group_by(Subject) %>% slice(1) %>% 
  select(-c("projectid","project","studyid","environmentName","subjectId","StudySiteId","SDVTier","siteid","Site","SiteNumber","SiteGroup","instanceId","InstanceName","InstanceRepeatNumber","folderid","Folder","FolderName","FolderSeq","TargetDays","DataPageId","PageRepeatNumber","RecordDate","RecordId","RecordPosition","MinCreated","MaxUpdated","SaveTS","StudyEnvSiteNumber")) %>%
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

label <- read_excel("Interim2_labels.xlsx", skip = 1)

label1 = label[!is.na(label[,3]),1:3]
colnames(label1)=c("Subject","Label","Batch")

label2 = label[!is.na(label[,6]),4:6]
colnames(label2)=c("Subject","Label","Batch")

label3 = label[!is.na(label[,9]),7:9]
colnames(label3)=c("Subject","Label","Batch")

Sepsis_Label <- rbind( rbind(label1 ,  label2), label3) 


res_t = data_join %>% rename("标本编号"="Patient Information_DXHSMP") %>%
       # unite(标本编号, c("StudyEnvSiteNumber", "ENROLL_NUM"),sep = "") %>% 
       #left_join(mdw_filter, by = "标本编号") %>%
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
      select(matches("Site")[1],Subject,标本编号,Time_Check,CBCADAT,Time,Enrollment_ENROLLYN_STD,Label,Batch,Label2,`Presenting Symptoms/Complaints (including symptom duration and intervention)_SYMOTH`,Lymph_index,Diff_MDW_Value,"CEC Adjudicator 1_SFDIAGA","CEC Adjudicator 1_FSDIAGARB","CEC Adjudicator 2_SFDIAGA","CEC Adjudicator 2_FSDIAGARB","CEC Arbitrator_SFDIAGA","CEC Arbitrator_FSDIAGARB") 

colnames(res)=c("Site","Subject", "标本编号","仪器分析时间检查","全血细胞分类计数分析日期时间","仪器真实分析时间","入组","剔除标签","是否更新入组策略","病毒感染提示","既有状况","Lymph_index","MDW", "Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3")

write.xlsx(arrange(res, Subject),  "Sepsis_STAT_new.xlsx",  colNames = TRUE)


######explore
Sepsis_type=c("脓毒性休克（严重脓毒症+低血压）","严重脓毒症（器官功能障碍或组织灌注不足）","脓毒症（全身炎症反应综合征+感染）")

res2 = res %>% filter(is.na(剔除标签)) %>%
       mutate_at(c('MDW'), as.numeric) %>% 
       drop_na(Adjudicator1_Sepsis2) %>% drop_na(Adjudicator2_Sepsis2) %>%
       mutate(Sepsis2=ifelse(Adjudicator1_Sepsis2 == Adjudicator2_Sepsis2, Adjudicator1_Sepsis2, Arbitrator_Sepsis2))%>% 
       drop_na(Sepsis2) %>% 
       mutate(Sepsis2_Sum = ifelse(Sepsis2 %in% Sepsis_type, "Sepsis" , Sepsis2))

my_comparisons <- list( c("感染，非脓毒症", "Sepsis"), c("感染，非脓毒症", "SIRS （≥2 SIRS标准）"), c("非全身炎症反应综合征/非感染（对照病例）","感染，非脓毒症") )


outpdf=paste("Sepsis","_cor.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=Lymph_index)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$Lymph_index), linetype = 2) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")


ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=MDW)) +
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

res3 =  res2 %>% filter(Sepsis2_Sum== "非全身炎症反应综合征/非感染（对照病例）") %>%
        mutate(MDW_Diag = ifelse(MDW > 20.5, "假阳", NA)) %>% 
        mutate(Ly_Diag = ifelse(Lymph_index > 11.68, "病毒感染", NA)) %>%
        select(MDW,MDW_Diag,Lymph_index,Ly_Diag,既有状况,Sepsis2_Sum,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 
        
write.xlsx(res3,  "Non-SIRS_screen.xlsx",  colNames = TRUE)


res4 =  res2 %>% filter(Sepsis2_Sum== "SIRS （≥2 SIRS标准）") %>%
        mutate(MDW_Diag = ifelse(MDW > 20.5, "假阳", NA)) %>% 
        mutate(Ly_Diag = ifelse(Lymph_index > 11.68, "病毒感染", NA)) %>%
        select(MDW,MDW_Diag,Lymph_index,Ly_Diag,既有状况,Sepsis2_Sum,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 
        
write.xlsx(res4,  "SIRS_screen.xlsx",  colNames = TRUE)