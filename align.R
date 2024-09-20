library(dplyr)
library(xlsx)
library(readxl)
library(compareDF)


two_a <- function(file){

mat1 = read_excel(file,sheet = 1,col_names = T)
mat2 = read_excel(file,sheet = 2,col_names = T)

inter_samples=intersect(mat1$筛选号,mat2$筛选号)
create_output_table(compare_df(arrange(mat1 %>% filter(筛选号 %in% inter_samples) , 筛选号),arrange(mat2 %>% filter(筛选号 %in% inter_samples) , 筛选号)), output_type = 'xlsx', file_name = paste0(file,"_queries.xlsx"))

print(inter_samples)

}

dir()[1]  %>% two_a()