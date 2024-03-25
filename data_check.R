library(dplyr)
library(xlsx)
library(readxl)


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