library(dplyr)
library(xlsx)
library(readxl)
library(compareDF)

dt = read_excel("CHN-189_DSF_CRA_20240320.xlsx",sheet = 2,col_names = T,skip=3)[-1,]

ins = read_csv("CSV-20240314.csv",skip=4)

tg_ins = ins %>% filter(检测名=="ThygIUO") %>% filter(grepl("^H",样本编号))  %>% select(样本编号,用户结果)


dt_query = dt %>% select(受试者入组号,`结果（ng/mL）`)%>% 
rename("样本编号"="受试者入组号") %>% 
full_join(tg_ins, by = "样本编号") %>%
mutate_at(c("结果（ng/mL）"), as.numeric)  %>%
mutate(用户结果=round(用户结果,2)) %>%
mutate(Query = ifelse(`结果（ng/mL）`!= 用户结果, "不等", "")) 


write_excel_csv(dt_query,  "TG_query.csv")







用户结果
table(dt$组别)


dt %>% filter(!grepl("^S{1}\\d{6}$", 筛选号))

dtt = dt %>% filter(!grepl("^P{1}\\d{3}$", 样本编号))
table(dtt$样本编号)

table(dt$性别)

sapply(dt$`年龄`, is.numeric)

dtt = dt %>% filter(!grepl("^\\d+", `结果（ng/L）`))
table(dtt$`结果（ng/L）`)

############
dt = read_excel(dir()[1],sheet = 2,col_names = T,skip=4)[-1,]

table(dt$性别)

sapply(dt$`年龄
（岁）`, is.numeric)

range(as.numeric(dt$对比试剂),na.rm =T)

range(as.numeric(dt$考核试剂),na.rm =T)



#####two group

ins = read_csv("DXI9000 Results_CEA and TPOAb_CSV_50%.csv")

predicate = read_csv("A2 Results_CEA and TPOAb_CSV_50%.csv",quote = "") %>% 
            rename(SampleID="Sample ID")

#######CEA
mat1 = ins %>% filter(TestName=="CEA") %>% 
            filter(SampleID %in% 219001:219050)  %>% 
            left_join(predicate,by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

mat2 = read_excel("CHN-219_CEA_DSF_CRA_TYH_20240418.xlsx", skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)


inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0("","_queries.xlsx"))



#######TPOAb
marker="TPOAb"

mat1 = ins %>% filter(TestName=="TPOAb") %>% 
            filter(SampleID %in% 228001:228050)  %>% 
            left_join(predicate  %>% filter(`Test Name`=="TPOAb"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

mat2 = read_excel("CHN-228_TPOAb_DSF_CRA_TYH_20240417.xlsx", skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（IU/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（IU/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))



#######DHEA-S

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")

marker=c("DHEA-S_IUO","DHEA-S")

mat1 = ins %>% filter(TestName%in%marker) %>% 
            filter(SampleID %in% c(170001:170061,170062:170136))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="DHE-S"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

mat2 = read_excel("CHN-170_DHEA-S_DSF_CRA_ZXY_20240828.xlsx", skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（µg/dL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（µg/dL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker[1],"_queries.xlsx"))


#######Total T3
marker="Total T3"
dsf="CHN-217_TT3_DSF_CRA_ZXY_20240531.xlsx"

mat1 = ins %>% filter(TestName==marker) %>% 
            filter(SampleID %in% 217001:217050)  %>% 
            left_join(predicate  %>% filter(`Test Name`=="TotT3"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))


#######hGH
marker="hGH"
dsf="CHN-227_hGH_DSF_CRA_ZXY_20240531.xlsx"

mat1 = ins %>% filter(TestName==marker) %>% 
            filter(SampleID %in% 227001:227051)  %>% 
            left_join(predicate  %>% filter(`Test Name`=="hGH2"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))



#######OV125Ag
marker="OV125Ag"
dsf="CHN-229_CA125_DSF_CRA_ZXY_20240529.xlsx"

mat1 = ins %>% filter(TestName==marker) %>% 
            filter(SampleID %in% 229001:229251)  %>% 
            left_join(predicate  %>% filter(`Test Name`==marker),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

write_excel_csv(ins %>% filter(TestName==marker) %>% select(SampleID,TestName,DoseResult),  paste0(marker,"_wrongnumbers.csv"))

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（U/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（U/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))



####batch3

#######CEA

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("CEA","CEA_IUO")
dsf="CHN-219_CEA_DSF_CRA_ZXY_20240627.xlsx"

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% 219001:219121)  %>% 
            left_join(predicate  %>% filter(`Test Name`=="CEA2"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

write_excel_csv(ins %>% filter(TestName==marker) %>% select(SampleID,TestName,DoseResult),  paste0(marker,"_wrongnumbers.csv"))

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))



#######BR15-3Ag

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("BR15-3Ag","BR15-3Ag_IUO")

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(230001:230046,228047:228069,230092:230103))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="BR15-3Ag"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName%in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))

mat2 = read_excel(dir()[2], skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（U/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（U/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))



#######TT3

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("Total T3","Total T3_IUO")
dsf="CHN-217_TT3_DSF_CRA_ZXY_20240809.xlsx"

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(217001:217056,217057:217114))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="TotT3"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

write_excel_csv(ins %>% filter(TestName==marker) %>% select(SampleID,TestName,DoseResult),  paste0(marker,"_wrongnumbers.csv"))

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))




#######CA125

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("OV125Ag","OV125Ag_IUO")
dsf="CHN-229_CA125_DSF_CRA_ZXY_20240823.xlsx"

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(229001:229070,229071:229119))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="OV125Ag"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

write_excel_csv(ins %>% filter(TestName==marker) %>% select(SampleID,TestName,DoseResult),  paste0(marker,"_wrongnumbers.csv"))

mat2 = read_excel(dsf, skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（U/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（U/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(marker,"_queries.xlsx"))



#######hGH

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("hGH","hGH_IUO")

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(227001:227058,227059:227117))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="hGH2"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))

mat2 = read_excel(dir()[2], skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))


#######Myoglobin

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("Myoglobin_IUO")

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(173001:173073))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="MYO"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))

mat2 = read_excel(dir()[2], skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)



inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))


#######Insulin

ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("Insulin_IUO")

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(226001:226079))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="Insulin"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))


mat2 = read_excel(dir()[2], skip = 2) %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ µIU/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ µIU/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))

#######Ferritin


ins = read_csv(dir()[3])

predicate = read_csv(dir()[1],quote = "") %>% 
            rename(SampleID="Sample ID")


markers=c("Ferritin_IUO")

mat1 = ins %>% filter(TestName %in% markers) %>% 
            filter(SampleID %in% c(159001:159063))  %>% 
            left_join(predicate  %>% filter(`Test Name`=="Ferritin"),by="SampleID") %>% 
            select(SampleID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))


mat2 = read_excel(dir()[2], skip = 2,col_types ="text") %>% 
          rename(SampleID="唯一可溯源编号") %>% 
          rename(DoseResult="DxI 9000仪器检测结果\r\n检测结果（ng/mL）") %>%
          rename(Result="Access 2仪器检测结果\r\n检测结果（ng/mL）") %>%
          select(SampleID,DoseResult,Result) %>%
          mutate_at(c('SampleID'), as.character)

inter_samples=intersect(mat1$SampleID,mat2$SampleID)

create_output_table(compare_df(arrange(mat1 %>% filter(SampleID %in% inter_samples) , SampleID),arrange(mat2 %>% filter(SampleID %in% inter_samples) , SampleID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))
