library(dplyr)
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(skimr)
library(data.table)

read_pre <- function(file,filter){
  print(file)
  edcform = read_csv(file,show_col_types = FALSE,guess_max = 100000 ) %>% select(-filter)
  
  if (edcform$DataPageName[1] == "实验室值：（在ED就诊后12小时内）"){
   edcform = edcform %>% select(-c("ANALYTE_STD")) %>% group_by(Subject) %>%
             pivot_wider(names_from = ANALYTE, values_from = c(LBDATE01,LBDATE02,LBDATE03,LBDATE04,LB01,LB02,LB03,LB04))

  }else if(edcform$DataPageName[1] == "Arterial blood gas (ABG): (within 12 hours from ED presentation)"){
   edcform = edcform %>% select(-c("ABGTERM_STD")) %>% group_by(Subject) %>%
             pivot_wider(names_from = ABGTERM, values_from = c(ABGDAT01,ABGDAT02,ABG01,ABG02))

  }else if(edcform$DataPageName[1] == "Vital signs: (within 12 hours from ED presentation)"){
   edcform = edcform %>% select(-c("VSTEST_STD")) %>% group_by(Subject) %>%
             pivot_wider(names_from = VSTEST, values_from = c(RSLTMAX,RSLTMIN))

  }else if(edcform$DataPageName[1] == "IV fluid administration and Use of Vasopressors (within 24 hours from ED presentation):"){
   edcform = dcast(setDT(edcform), Subject + DataPageName ~ paste0("IV_", rowid(Subject)), value.var = c('IV_VOL','IV_DAT')) %>% 
             full_join(dcast(setDT(edcform), Subject ~ paste0("VAS_", rowid(Subject)), value.var = c('VAS_DUR','VAS_DAT')) ,by = "Subject")

  }else if(edcform$DataPageName[1] == "Bacterial Testing Log"){
    source_term = setdiff(colnames(edcform),c("Subject","DataPageName","SPECIMYN","SPECIMYN_STD","SEROYN","SEROYN_STD","SERORSLT","SERORSLT_STD","SEROULAP","SEROUSAP","SEROSPC","SEROCUL","SEROCUL_STD","SOURCE","SOURCE_STD")) 
    edcform = edcform %>% group_by(Subject) %>%
             pivot_wider(names_from = SOURCE, values_from = source_term) %>% 
             distinct(Subject, .keep_all=TRUE)
  }
  else{print(edcform %>% group_by(Subject) %>% filter(n() > 1))}
  edcform %>% rename_with(~ paste0(edcform$DataPageName[1],"_", .), -c("Subject","DataPageName"))
}


####SAS On Demand
long_filter_terms=c("projectid","project","studyid","environmentName","subjectId","StudySiteId","SDVTier","siteid","Site","SiteNumber","SiteGroup","instanceId","InstanceName","InstanceRepeatNumber","folderid","Folder","FolderName","FolderSeq","TargetDays","DataPageId","PageRepeatNumber","RecordDate","RecordId","RecordPosition","MinCreated","MaxUpdated","SaveTS","StudyEnvSiteNumber")

data_join <- list.files(path = "./CHN_097_EDC/", # Identify all CSV files
                       pattern = "*.CSV|*.txt", full.names = TRUE) %>% 
  lapply(read_pre,filter=long_filter_terms) %>%                              # Store all files in list
  reduce(full_join, by = "Subject")                      # Full-join data sets into one data set 


####data listing
short_filter_terms = c("project","studyid","SDVTier","SiteNumber","SiteGroup","InstanceName","InstanceRepeatNumber","FolderName","RecordPosition","SaveTS","StudyEnvSiteNumber")
lising_join = list.files(path = "./listing/", # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_pre,filter=short_filter_terms) %>%                              # Store all files in list
  reduce(full_join, by = "Subject")


###instrument file
ins <- read_csv("Interim2_instruments.csv")

instrument_filter = ins %>% 
select(RDFilename,Wbc,Ne_pct,Mo_pct,matches("MDW"),Ly_no,Ly_pct,matches("Ly_DC"),matches("Ly_Op_Mean")) %>% 
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
       left_join(instrument_filter, by = "标本编号") %>% 
       left_join(Sepsis_Label, by = "Subject") 


###generate final table
# res_t$Duplicate <- duplicated(res_t$Subject)
res_dat = res_t %>% drop_na(Subject) %>% 
      rename(CBCADAT="Blood Collection for the Study_CBCADAT") %>%
      mutate(CBCADAT=format(strptime(CBCADAT, "%d/%b/%Y %H:%M"),format='%Y-%m-%d %H:%M')) %>% 
      mutate(Time=format(as_datetime(Time),format='%Y-%m-%d %H:%M')) %>% 
      mutate(Time_Check = ifelse(CBCADAT==Time, "Yes", "No")) %>%
      group_by(Subject) %>% mutate(freq=n()) %>%
      filter(!(freq!=1 & CBCADAT!=Time)) %>%
      mutate(Label = ifelse(Enrollment_ENROLLYN_STD=="N" & is.na(Label), "未入组", Label))


####for jimson

trans_chi <- function(x){
    case_when(x=="非全身炎症反应综合征/非感染（对照病例）" ~ "Non-SIRS/non-infection (control case)",
              x=="SIRS （≥2 SIRS标准）"   ~ "SIRS (>= 2 SIRS criteria)" ,
              x=="感染，非脓毒症"   ~ "Infection without sepsis",
              x=="脓毒症（全身炎症反应综合征+感染）"    ~ "Sepsis (2 SIRS + infection)",
              x=="脓毒性休克（严重脓毒症+低血压）"    ~ "Septic shock (severe sepsis + hypotension)" ,
              x=="严重脓毒症（器官功能障碍或组织灌注不足）"   ~ "Severe sepsis (organ dysfunction or tissue hypoperfusion)",
              x=="非脓毒症（病例对照）" ~ "Non-sepsis (case control)",
              x=="确认感染，无脓毒症" ~ "Confirmed infection without sepsis",
              x=="脓毒症（SOFA 增加>=2）" ~ "Sepsis (SOFA increase >= 2)",
              x=="脓毒性休克（脓毒症+低血压）" ~ "Septic shock (sepsis + hypotension)"
              )
}

trans_site <- function(x){
    case_when(x=="01" ~ "Peking Union Medical College Hospital"  ,
              x=="02"   ~ "West China School of Medicine; Sichuan University" ,
              x=="03"   ~ "Zhejiang University College of Medicine- Second Affiliated Hospital",
              )
}

res_stat = res_dat %>% 
      select(matches("Site")[1],Subject,"Demographic Information of Patients_SEX_STD","Demographic Information of Patients_AGE",Label,Batch,Diff_MDW_Value,"CEC Adjudicator 1_SFDIAGA","CEC Adjudicator 1_FSDIAGARB","CEC Adjudicator 2_SFDIAGA","CEC Adjudicator 2_FSDIAGARB","CEC Arbitrator_SFDIAGA","CEC Arbitrator_FSDIAGARB") 

colnames(res_stat)=c("Site","Subject", "性别","年龄", "剔除标签","是否更新入组策略","MDW", "Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3")

res_stat2 = res_stat %>% mutate_at("Site",trans_site) %>% 
            mutate_at(c("Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3"),trans_chi) 

write.xlsx(arrange(res_stat2, Subject),  "Sepsis_STAT.xlsx",  colNames = TRUE)

##check
#res_com = res_stat2 %>%  right_join(lising_join, by = "Subject") %>% 
#          mutate(Adjudicator1_Sepsis2_check = ifelse(Adjudicator1_Sepsis2 != `CEC Adjudicator 1_SFDIAGA`, "Wrong" , "")) %>%
#          mutate(Adjudicator1_Sepsis3_check = ifelse(Adjudicator1_Sepsis3 != `CEC Adjudicator 1_FSDIAGARB`, "Wrong" , "")) %>%
#          select("Site","Subject","Adjudicator1_Sepsis2_check","Adjudicator1_Sepsis2","CEC Adjudicator 1_SFDIAGA","Adjudicator1_Sepsis3_check","Adjudicator1_Sepsis3","CEC Adjudicator 1_FSDIAGARB")
#write.xlsx(arrange(res_com, Subject),  "Sepsis_com_test.xlsx",  colNames = TRUE)

#report version

 
res = res_dat  %>% 
     #compute Lymph_index
     mutate_at(c('Ly_DC_Mean', 'Ly_DC_SD', 'Ly_Op_Mean'), as.numeric) %>% 
     mutate(Lymph_index =Ly_DC_Mean*Ly_DC_SD/Ly_Op_Mean) %>% 
     mutate(Ly_Label = ifelse(Lymph_index > 11.68, "Virus_Infection_Lymph_index", NA)) %>% 
     #仪器血常规
     mutate_at(c('Wbc','Ne_pct','Ly_no','Ly_pct'), as.numeric) %>% 
     mutate(BV_Label = ifelse(Ly_pct > 50, "Virus_Infection_IN_Blood_routine", NA)) %>%
     mutate(BB_Label = ifelse(Wbc > 9.5 & Ne_pct > 70, "Bac_Infection_IN_Blood_routine", NA)) %>%
     #EDC血常规
     rename(EDC_WBC1="实验室值：（在ED就诊后12小时内）_LB01_白细胞（×10^9/L)")%>%
     rename(EDC_Ly_pct1="实验室值：（在ED就诊后12小时内）_LB01_淋巴细胞 %")%>%
     rename(EDC_Ne_pct1="实验室值：（在ED就诊后12小时内）_LB01_中性粒细胞 %")%>%
     mutate(across(c("EDC_WBC1","EDC_Ly_pct1","EDC_Ne_pct1"), readr::parse_number))%>%
     mutate(EDC_Ly_num1 = EDC_WBC1*EDC_Ly_pct1)%>%
     #EDC细菌指标
     rename(EDC_IL6_1="实验室值：（在ED就诊后12小时内）_LB01_白细胞介素6 (pg/mL) (如适用)")%>%
     rename(EDC_PCT1="实验室值：（在ED就诊后12小时内）_LB01_降钙素原 （ng/mL）")%>%
     rename(EDC_CRP1="实验室值：（在ED就诊后12小时内）_LB01_C反应蛋白 （mg/L）")%>%
     rename(EDC_INR1="实验室值：（在ED就诊后12小时内）_LB01_国际标准化比值 (INR)") %>%
     mutate(across(c("EDC_PCT1","EDC_CRP1","EDC_IL6_1","EDC_INR1"), readr::parse_number))%>%
     mutate(IL6_Label = ifelse(EDC_IL6_1 > 7, "Bac_Infection_IL6", NA)) %>% 
     mutate(PCT_Label = ifelse(EDC_PCT1 > 1.1, "Bac_Infection_PCT", NA)) %>% 
     mutate(CRP_Label = ifelse(EDC_CRP1 > 100, "Bac_Infection_CRP", NA)) %>%
     #LIP score
     mutate(LIP_L = cut(Ly_no, breaks = c(0, 0.7, 1), right = FALSE, labels = c("2", "1"))) %>%
     mutate(LIP_I = cut(EDC_INR1, breaks = c(1.2,1.4,Inf), right = FALSE, labels = c("1", "2"))) %>%
     mutate(LIP_P = cut(EDC_PCT1, breaks = c(0.5,2,Inf), right = FALSE, labels = c("1", "2"))) %>%
     mutate_at(c('LIP_L', 'LIP_I', 'LIP_P'), as.character) %>%
     mutate_at(c('LIP_L', 'LIP_I', 'LIP_P'), as.numeric) %>% 
     mutate(LIP = LIP_L + LIP_I + LIP_P) %>%
     mutate(LIP_Label = ifelse(LIP > 2, "Sepsis3", NA)) %>%
     # test wrong
     # filter(Subject=="01-0244") %>%
     # select(Subject,LIP,LIP_L,Ly_no,LIP_I,EDC_INR1 ,LIP_P,EDC_PCT1)
     #EDC information
     #既有状况_感染影响
     rename(IF_Ac="Pre-existing Conditions_ED08OTH") %>%
     rename(IF_Ac_YN="Pre-existing Conditions_ED08_STD") %>%
     rename(IF_Pre="Pre-existing Conditions_ED09OTH") %>%
     rename(IF_Pre_YN="Pre-existing Conditions_ED09_STD") %>%
     #既有状况_免疫抑制
     rename(ImSu_Chemo="Pre-existing Conditions_ED04OTH") %>%
     rename(ImSu_Chemo_YN="Pre-existing Conditions_ED04_STD") %>%
     rename(ImSu_Long="Pre-existing Conditions_ED05OTH") %>%
     rename(ImSu_Long_YN="Pre-existing Conditions_ED05_STD") %>%
     rename(ImSu_Other_YN="Pre-existing Conditions_ED03_STD") %>%
     rename(ImSu_HU_YN="Pre-existing Conditions_ED06_STD") %>%
     #既有状况_血细胞异常
     rename(BC_Ab="Pre-existing Conditions_ED07OTH") %>%
     rename(BC_Ab_YN="Pre-existing Conditions_ED07_STD") %>%
     rename(BC_Ne_De_YN="Pre-existing Conditions_ED01_STD") %>%
     rename(BC_Ne_Th_YN="Pre-existing Conditions_ED02_STD") %>%
     mutate(ED_IF = case_when(IF_Ac_YN=="Y" ~ "Infection",
                              IF_Pre_YN=="Y" ~ "Infection")) %>%
     mutate(ED_ImSu = case_when(ImSu_Chemo_YN=="Y" ~ "Immu_Supression",
                                ImSu_Long_YN=="Y" ~ "Immu_Supression",
                                ImSu_Other_YN=="Y" ~ "Immu_Supression",
                                ImSu_HU_YN=="Y" ~ "Immu_Supression")) %>%  
    mutate(ED_BC = case_when(BC_Ab_YN=="Y" ~ "Blood_Cell_Abnormal",
                             BC_Ne_De_YN=="Y" ~ "Blood_Cell_Abnormal",
                             BC_Ne_Th_YN=="Y" ~ "Blood_Cell_Abnormal")) %>%                           
     #既往病史
     rename(PM_Gastr="Past Medical History_GASTRYN_STD") %>%
     rename(PM_Hema="Past Medical History_HEMATYN_STD") %>%
     rename(PM_AutoIm="Past Medical History_AUTOIMYN_STD") %>%
     rename(PM_Car="Past Medical History_CARDIYN_STD") %>%
     rename(PM_Genit="Past Medical History_GENITYN_STD") %>%
     rename(PM_Respi="Past Medical History_RESPIYN_STD") %>%
     rename(PM_Meta="Past Medical History_METAYN_STD") %>%
     rename(PM_Cns="Past Medical History_CNSYN_STD") %>%
     rename(PM_Renal="Past Medical History_RENALYN_STD") %>%
     rename(PM_Hepati="Past Medical History_HEPATIYN_STD") %>%
     rename(PM_Other="Past Medical History_ORTOHSPC") %>%
     #感染处理
     rename(IW="Infection Workup_IW_YN") %>%
     rename(SEROYN="Bacterial Testing Log_SEROYN") %>%
     rename(VIRALRES="Other Culture and Rapid Test_VIRALRES") %>%
     rename(FUNGALRE="Other Culture and Rapid Test_FUNGALRE") %>%
     rename(OTHTPERF="Other Culture and Rapid Test_OTHTPERF") %>%
     # 抗生素抗病毒药物
     rename(ANTIBYN="Antibiotics, antiviral and antifungal drugs_ANTIBYN") %>%
     rename(ANTIVYN="Antibiotics, antiviral and antifungal drugs_ANTIVYN") %>%
     rename(ANTIFYN="Antibiotics, antiviral and antifungal drugs_ANTIFYN") %>%
     rename(OTHINYN="Antibiotics, antiviral and antifungal drugs_OTHINYN") %>%
     # 诊断手术
     rename(TOMO_INF="Surgery/Diagnostics_TOMO_INF") %>%
     rename(TOMOFIND="Surgery/Diagnostics_TOMOFIND") %>%
     rename(ABCTFIND="Surgery/Diagnostics_ABCTFIND") %>%
     rename(ABULFIND="Surgery/Diagnostics_ABULFIND") %>%
     rename(OTHEXFIN="Surgery/Diagnostics_OTHEXFIN") %>%
     rename(SURGPFIN="Surgery/Diagnostics_SURGPFIN") %>%
     rename(XRAYFIND="Surgery/Diagnostics_XRAYFIND") %>%
     select(matches("Site")[1],Subject,标本编号,Time_Check,CBCADAT,Time,Enrollment_ENROLLYN_STD,Label,Batch,IW,SEROYN,ANTIBYN,VIRALRES,ANTIVYN,FUNGALRE,ANTIFYN,OTHTPERF,OTHINYN,IL6_Label,EDC_IL6_1,PCT_Label,EDC_PCT1,CRP_Label,EDC_CRP1,BV_Label,Wbc,EDC_WBC1,Ly_pct,EDC_Ly_pct1,BB_Label,Ne_pct,EDC_Ne_pct1,Ly_Label,Lymph_index,Diff_MDW_Value,LIP_Label,LIP,Ly_no,EDC_INR1,ED_IF,IF_Ac,IF_Pre,ED_ImSu,ImSu_Chemo,ImSu_Long,ImSu_Other_YN,ImSu_HU_YN,ED_BC,BC_Ab,BC_Ne_De_YN,BC_Ne_Th_YN,"CEC Adjudicator 1_SFDIAGA","CEC Adjudicator 1_FSDIAGARB","CEC Adjudicator 2_SFDIAGA","CEC Adjudicator 2_FSDIAGARB","CEC Arbitrator_SFDIAGA","CEC Arbitrator_FSDIAGARB",`Presenting Symptoms/Complaints (including symptom duration and intervention)_SYMOTH`,PM_Gastr,PM_Hema,PM_AutoIm,PM_Car,PM_Genit,PM_Respi,PM_Meta,PM_Cns,PM_Renal,PM_Hepati,PM_Other,TOMO_INF,TOMOFIND,ABCTFIND,ABULFIND,OTHEXFIN,SURGPFIN,XRAYFIND) 
 
 res_c=res
 colnames(res_c)=c("Site","Subject", "标本编号","仪器分析时间检查","全血细胞分类计数分析日期时间","仪器真实分析时间","入组","剔除标签","是否更新入组策略","感染处理","细菌检测","抗生素","病毒检测","抗病毒","真菌检测","抗真菌","其他感染","其他抗感染药物", "EDC_IL6_Bacteria_Infection","EDC_IL6_1","EDC_PCT_Bacteria_Infection","EDC_PCT1","EDC_CRP_Bacteria_Infection","EDC_CRP1","血常规_病毒感染提示","INS_WBC","EDC_WBC1","INS_Ly_Percent","EDC_Ly_Percent1","血常规_细菌感染提示","INS_Ne_Percent","EDC_Ne_Percent1","淋巴指数_病毒感染提示","Lymph_index","MDW","LIP_Sepsis3","LIP_score","Ly_Num","EDC_INR1","既有状况_感染","活动性感染而接受抗生素治疗","抗生素初级预防治疗","既有状况_免疫抑制","化疗导致免疫抑制","免疫抑制剂长期治疗","因HIV or 器官或骨髓移植而疑似免疫抑制","羟基脲治","既有状况_血细胞异常","患有影响全血细胞分类计数的血液病受试者","中性粒细胞减少症","受试者接受中性粒细胞减少症治疗","Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3","出现症状" ,"PM_Gastr","PM_Hema","PM_AutoIm","PM_Car","PM_Genit","PM_Respi","PM_Meta","PM_Cns","PM_Renal","PM_Hepati","PM_Other","计算机断层扫描浸润或固结","计算机断层扫描检查结果","腹部计算机断层扫描检查结果","腹部超声检查结果","其他相关检查检查结果","手术检查结果","胸部x光检结果")
 
 write.xlsx(arrange(res_c, Subject),  "Sepsis_CRA3.xlsx",  colNames = TRUE)


######explore

colnames(res)=c("Site","Subject", "标本编号","仪器分析时间检查","全血细胞分类计数分析日期时间","仪器真实分析时间","入组","剔除标签","是否更新入组策略","感染处理","细菌检测","抗生素","病毒检测","抗病毒","真菌检测","抗真菌","其他感染","其他抗感染药物","EDC_IL6_Bacteria_Infection","EDC_IL6_1","EDC_PCT_Bacteria_Infection","EDC_PCT1","EDC_CRP_Bacteria_Infection","EDC_CRP1","BV_Label","INS_WBC","EDC_WBC1","INS_Ly_Percent","EDC_Ly_Percent1","BB_Label","INS_Ne_Percent","EDC_Ne_Percent1","Ly_Label","Lymph_index","MDW", "LIP_Sepsis3","LIP_score","Ly_Num","EDC_INR1","ED_IF","IF_Ac","IF_Pre","ED_ImSu","ImSu_Chemo","ImSu_Long","ImSu_Other_YN","ImSu_HU_YN","ED_BC","BC_Ab","BC_Ne_De_YN","BC_Ne_Th_YN","Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3","出现症状" ,"PM_Gastr","PM_Hema","PM_AutoIm","PM_Car","PM_Genit","PM_Respi","PM_Meta","PM_Cns","PM_Renal","PM_Hepati","PM_Other","计算机断层扫描浸润或固结","计算机断层扫描检查结果","腹部计算机断层扫描检查结果","腹部超声检查结果","其他相关检查检查结果","手术检查结果","胸部x光检结果")
 
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

ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=MDW)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$MDW), linetype = 2) +
stat_compare_means( ref.group = ".all.",method = "wilcox.test", label = "p.signif")

ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=MDW)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$MDW), linetype = 2) +
stat_compare_means( ref.group = "非全身炎症反应综合征/非感染（对照病例）",method = "wilcox.test", label = "p.signif")


ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=Lymph_index)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$Lymph_index), linetype = 2) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")


ggplot(data=res2,aes(x=factor(Sepsis2_Sum,level = c("非全身炎症反应综合征/非感染（对照病例）","SIRS （≥2 SIRS标准）", "感染，非脓毒症", "Sepsis")),y=MDW)) +
geom_boxplot()+
geom_jitter(width = 0.2, alpha = 0.5, color = 'red') +
geom_hline(yintercept = mean(res2$MDW), linetype = 2) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")


ggplot(res2, aes(x = MDW, y = Lymph_index)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "blue", alpha = 0.3) +
  theme_bw()+stat_cor(data=res2, method = "spearman")

ggplot(res2, aes(x = MDW, y = Lymph_index, color=Sepsis2_Sum)) + geom_point() + stat_smooth(method = 'lm')

dev.off()

res3 =  res2 %>% filter(Sepsis2_Sum== "非全身炎症反应综合征/非感染（对照病例）") %>%
        mutate(MDW_Diag = ifelse(MDW > 20, "MDW>20", NA)) %>% 
        select(Sepsis2_Sum,,ED_IF,IF_Ac,IF_Pre,ED_ImSu,ImSu_Chemo,ImSu_Long,ImSu_Other_YN,ImSu_HU_YN,ED_BC,BC_Ab,BC_Ne_De_YN,BC_Ne_Th_YN,MDW,MDW_Diag,LIP_Sepsis3,LIP_score,"感染处理","细菌检测","抗生素","病毒检测","抗病毒","其他感染",Ly_Label,BV_Label,EDC_IL6_Bacteria_Infection,EDC_PCT_Bacteria_Infection,EDC_CRP_Bacteria_Infection,BB_Label,出现症状,"其他相关检查检查结果","手术检查结果",PM_Gastr,PM_Hema,PM_AutoIm,PM_Car,PM_Genit,PM_Respi,PM_Meta,PM_Cns,PM_Renal,PM_Hepati,PM_Other,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 


write.xlsx(res3,  "Non-SIRS_screen.xlsx",  colNames = TRUE)


res4 =  res2 %>% filter(Sepsis2_Sum== "SIRS （≥2 SIRS标准）") %>%
        mutate(MDW_Diag = ifelse(MDW > 20, "MDW>20", NA)) %>% 
        select(Sepsis2_Sum,,ED_IF,IF_Ac,IF_Pre,ED_ImSu,ImSu_Chemo,ImSu_Long,ImSu_Other_YN,ImSu_HU_YN,ED_BC,BC_Ab,BC_Ne_De_YN,BC_Ne_Th_YN,MDW,MDW_Diag,LIP_Sepsis3,LIP_score,"感染处理","细菌检测","抗生素","病毒检测","抗病毒","其他感染",Ly_Label,BV_Label,EDC_IL6_Bacteria_Infection,EDC_PCT_Bacteria_Infection,EDC_CRP_Bacteria_Infection,BB_Label,出现症状,"其他相关检查检查结果","手术检查结果",PM_Gastr,PM_Hema,PM_AutoIm,PM_Car,PM_Genit,PM_Respi,PM_Meta,PM_Cns,PM_Renal,PM_Hepati,PM_Other,"Adjudicator1_Sepsis2", "Adjudicator1_Sepsis3","Adjudicator2_Sepsis2", "Adjudicator2_Sepsis3","Arbitrator_Sepsis2", "Arbitrator_Sepsis3") 
      
write.xlsx(res4,  "SIRS_screen.xlsx",  colNames = TRUE)


write.xlsx(rbind(res3,res4) %>% filter(!is.na(MDW_Diag)),  "Negative_screen2.xlsx",  colNames = TRUE)


print(paste0("非全身炎症反应综合征/非感染（对照病例）","Lymph_index Mean:", mean(res3$Lymph_index) %>% round(2) , "+_", sd(res3$Lymph_index) %>% round(2) ))     
print(paste0("SIRS （≥2 SIRS标准）","Lymph_index Mean:", mean(res4$Lymph_index)%>% round(2), "+-", sd(res4$Lymph_index)%>% round(2)))     


print(paste0("非全身炎症反应综合征/非感染（对照病例）","MDW:", mean(res3$MDW)%>% round(2), "+-", sd(res3$MDW)%>% round(2)))     
print(paste0("SIRS （≥2 SIRS标准）","MDW:", mean(res4$MDW)%>% round(2), "+-", sd(res4$MDW)%>% round(2)))     
