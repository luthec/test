library(readxl)
library(tableone)
library(dplyr)


empyrosis_data1 <- read_excel("burn/empyrosis_data1.xlsx", 
     col_types = c("text", "text", "numeric", 
         "numeric", "numeric", "text", "numeric", 
         "numeric", "numeric", "text", "numeric", 
         "numeric", "text", "numeric", "text", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "text", "text", "text", 
         "text", "numeric", "text", "numeric", 
         "numeric", "text", "text", "numeric", 
         "numeric"))


empyrosis_t =  empyrosis_data1 %>% select(c("体重（kg","BMI","TBSA","烧伤指数","@90天死亡（成活0死亡1）","呕吐（无0有1）"))

## List numerically coded categorical variables
#factorVars <- c("@90天死亡（成活0死亡1）","呕吐（无0有1）")
## Create a variable list. Use dput(names(pbc))
#vars <- c("体重（kg","BMI","TBSA","烧伤指数")

#tableOne <- CreateTableOne(data = empyrosis_t , strata = "呕吐（无0有1）", vars = vars, factorVars = factorVars, smd = TRUE)


tableOne <- CreateTableOne(vars = colnames(select(empyrosis_t, -"呕吐（无0有1）")), 
                           strata = c("呕吐（无0有1）"), 
                           data = empyrosis_t)

print(
  tableOne,
  nonnormal = c("TBSA","烧伤指数"),exact = c("@90天死亡（成活0死亡1）"),
  showAllLevels = TRUE)                           

ffff