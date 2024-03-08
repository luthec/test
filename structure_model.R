library(performance)
library(psych)

library(lavaan)
library(plspm)
library(semPlot)

library(tidyverse)
library(readxl)
library(dplyr)

test_datasets <- read_excel("test_datasets.xlsx", 
    col_types = c("numeric", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "text", "numeric", "numeric", "numeric", 
        "text", "text", "numeric", "text", 
        "numeric", "text", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "text", "numeric", 
        "text", "text", "text", "text", "text", 
        "numeric", "text", "numeric", "numeric", 
        "text", "text", "numeric", "numeric", 
        "numeric", "text", "text", "text", 
        "numeric", "text", "text", "text"))

dt_select <- test_datasets %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(HCTdALB = test_datasets$`@1d_hct`/test_datasets$`@1d_ALB`) %>%
   mutate(`@1d_ALB` = cut(`@1d_ALB`, breaks = c(-Inf, 25, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`24h血糖（mmol/L）` = cut(`24h血糖（mmol/L）`, breaks = c(-Inf, 14, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`= ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==0, "0", "1")) %>% rename("基础疾病"="基础疾病(0无1高血压2糖尿病3两者都有4其他)") %>% 
   mutate(`并发症`= ifelse(`并发症`==0, "0", "1")) %>% 
   select(-c("病案号","姓名","体重（kg","身高m","天数","中性粒细胞","吸入性损伤...46","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）"))   %>% 
   rename_with(~str_remove(., '[@）]+')) %>% 
   mutate_if(is.character, as.factor) %>%  
   rename("耐受"="主要结局：耐受0，不耐受1") %>% rename("吸入性损伤"="吸入性损伤...12") 

colnames(dt_select) <- gsub('[(（].*','',colnames(dt_select))

dt_train = dt_select  %>% 
    select(c("年龄","性别","BMI","基础疾病","TBSA","烧伤指数","III度","休克","吸入性损伤","鼻饲","启动时间","首次进水时间h","HCTdALB","1d_ALB","1dHB","1d_plt","1d_淋巴细胞","1d_TP","1d_TBIL","1d_DBIL","1d_CRE","1d_BUN","1d_PA","1d_LAC","耐受","呕吐","胃肠减压","腹胀","腹泻","24h血糖","首次排便时间","CRRT","并发症","脓毒症","90天死亡","气切","呼吸机","天数"))

colnames(dt_train)=c("Age_H","Sex","Bmi","Basic_Disease","TBSA","Burn_index","III_index","Shock","Inhalation_injury","Nasal","Start_Time_D","First_DrinkTime_H","HCTdALB","ALB","HB","Plt","Lymphocyte","TP","TBIL","DBIL","CRE","BUN","PA","LAC","Res","Emesis","Gastro_De","Abdominal_distension","Diarrhea","Glu_24h","F_defecation_D","CRRT","Complication","Sepsis","Dease","Gas_cut","Ventilator","Ventilator_Day")



LVdata.model <-'
# CFAs
Demo_index =~ Age_H + Sex + Bmi + Basic_Disease
Blood_index =~ HCTdALB + ALB + HB + Plt + Lymphocyte + TP + TBIL + DBIL + CRE + BUN + PA + LAC
Burn_indices =~ TBSA + Burn_index + III_index
Digestive_tract =~ Shock + Inhalation_injury + Nasal + Start_Time_D + First_DrinkTime_H
Prognosis =~ Emesis + Gastro_De + Abdominal_distension + Diarrhea + Glu_24h + F_defecation_D + CRRT + Complication + Sepsis + Dease
# Regressions
Res ~ Demo_index + Blood_index + Burn_indices + Digestive_tract
Prognosis ~ Res + Demo_index + Blood_index + Burn_indices + Digestive_tract
'
# Residuals
#y1 ~~ y5
#y2 ~~ y4 + y6
#y3 ~~ y7
#y4 ~~ y8
#y6 ~~ y8
#'

LVdata.model <-'
# CFAs
Demo_index =~ Age_H + Sex + Bmi + Basic_Disease
Blood_index =~ HCTdALB + ALB + HB + Plt + Lymphocyte + TP + TBIL + DBIL + CRE + BUN + PA + LAC
Burn_indices =~ TBSA + Burn_index + III_index 
Digestive_tract =~ Shock + Inhalation_injury + Nasal + Start_Time_D + First_DrinkTime_H
# Regressions
Res ~ Demo_index + Blood_index + Burn_indices + Digestive_tract
'
fit <- cfa(LVdata.model, data= dt_train, std.lv=T,ordered = c("Age_H", "Sex", "Basic_Disease", "ALB", "Shock", "Inhalation_injury", "Nasal", "Gastro_De", "Abdominal_distension", "Diarrhea", "Glu_24h", "CRRT", "Complication", "Sepsis", "Dease", "Res"))
summary(fit, fit.measures=TRUE)

LVdata.fit <- lavaan:::sem(LVdata.model, data= dt_train, std.lv=T,ordered = c("Age_H", "Sex", "Basic_Disease", "ALB", "Shock", "Inhalation_injury", "Nasal", "Gastro_De", "Abdominal_distension", "Diarrhea", "Glu_24h", "CRRT", "Complication", "Sepsis", "Dease", "Res"))

summary(LVdata.fit, fit.measures=TRUE, standardized = TRUE)


semPaths(LVdata.fit, "std", edge.label.cex = 0.5, curvePivot = TRUE,intercepts = FALSE)

semPaths(LVdata.fit, as.expression = c("nodes", "edges"), sizeMan = 3, sizeInt = 1, 
    sizeLat = 6, label.prop = 0.5, curve = 0.5,  groups = "latents", 
    intercepts = FALSE, borders = FALSE, label.norm = "O")