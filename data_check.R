library(dplyr)
library(xlsx)
library(readxl)


dt = read_excel(dir()[2],sheet = 2,col_names = T)

ins = read.csv(dir()[2], header = TRUE,skip=4)

tg_ins = ins %>% as.tibble() %>% filter(检测名=="ThygIUO") %>% filter(grepl("^H",样本编号))  %>% select(样本编号,用户结果)

用户结果
table(dt$组别)


dt %>% filter(!grepl("^S{1}\\d{6}$", 筛选号))

dtt = dt %>% filter(!grepl("^P{1}\\d{3}$", 样本编号))
table(dtt$样本编号)

table(dt$性别)

sapply(dt$`年龄`, is.numeric)

dtt = dt %>% filter(!grepl("^\\d+", `NT-proBNP 检测结果（ng/L）`))
table(dtt$`NT-proBNP 检测结果（ng/L）`)

############
dt = read_excel(dir()[1],sheet = 2,col_names = T,skip=4)[-1,]

table(dt$性别)

sapply(dt$`年龄
（岁）`, is.numeric)

range(as.numeric(dt$对比试剂),na.rm =T)

range(as.numeric(dt$考核试剂),na.rm =T)