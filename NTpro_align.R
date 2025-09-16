library(dplyr)
library(compareDF)
library(readxl)
library(tidyverse)
library(eatTools)

trans_name <- function(x){
    case_when(x=="非急性心衰" ~ "非AHF"  ,
              x=="急性心衰" ~ "AHF" ,
    )
}

Site1 <- read_excel("Site 1-A.xlsx") %>%
    # rename(贝克曼上机数据="贝克曼上机结果") %>% 
    mutate_at(c("贝克曼上机数据","心超LVEF","体重","eGFR","左室缩短分数","左室舒张末期内径","左室收缩期末期内径")  , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 1))) %>% suppressWarnings() %>%
    # rename(诊断结果="第三方专家确诊诊断") %>% 
    # mutate_at("诊断结果",trans_name)%>% 
    select(!matches("时间"))


Site2 <- read_excel("Site 2-B.xls") %>%
    rename(贝克曼上机数据="贝克曼上机结果") %>% 
    mutate_at(c("贝克曼上机数据","心超LVEF","体重","eGFR","左室缩短分数","左室舒张末期内径","左室收缩期末期内径")   , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 1))) %>% suppressWarnings() %>%
    # mutate_at("诊断结果",trans_name)%>% 
    select(!matches("时间"))



Site3 <- read_excel("Site3-C.xlsx") %>%
    rename(贝克曼上机数据="贝克曼上机结果") %>% 
    mutate_at(c("贝克曼上机数据","心超LVEF","体重","eGFR","左室缩短分数","左室舒张末期内径","左室收缩期末期内径")   , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 1))) %>% suppressWarnings() %>%
    # mutate_at("诊断结果",trans_name)%>% 
    select(!matches("时间"))


allsite = rbind_common(rbind_common(Site1, Site2),Site3)



database <- read_excel("SierraNT-proBNP_database.xlsx")%>% 
    select(!matches("入组")) %>% 
    mutate_at(c("贝克曼上机数据","心超LVEF","体重","eGFR","左室缩短分数","左室舒张末期内径","左室收缩期末期内径") , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 1))) %>% suppressWarnings() 

target_columns = intersect(database  %>% colnames(), allsite %>% colnames() )

final = database %>% select(target_columns) 
original = allsite  %>% select(target_columns) 

inter_samples=intersect(final$筛选号,original$筛选号)
compare_df(arrange(final, 筛选号), arrange(original %>% filter(筛选号 %in% inter_samples) ,筛选号)) %>%  create_output_table(output_type = 'xlsx', file_name = paste0("NT-PRO","_queries.xlsx"))

