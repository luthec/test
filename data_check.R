library(dplyr)
library(xlsx)
library(readxl)


dt = read_excel(dir()[1],sheet = 2,col_names = T)



table(dt$组别)


dt %>% filter(!grepl("^S{1}\\d{6}$", 筛选号))

dtt = dt %>% filter(!grepl("^P{1}\\d{3}$", 样本编号))
table(dtt$样本编号)

table(dt$性别)

sapply(dt$`年龄（岁）`, is.numeric)

dtt = dt %>% filter(!grepl("^\\d+", `NT-proBNP 检测结果（ng/L）`))
table(dtt$`NT-proBNP 检测结果（ng/L）`)