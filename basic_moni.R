library(readxl)
library(tidyverse)

table(CHN_184_7_test_results[,6])

uE3检测结果
    其他 检测失败 
     988       11 


table(CHN_184_7_test_results[,12])

最终检测结果
其他 
 999 

CHN_184_IE = CHN_184_4_6_PD_table %>% filter(是否符合剔除标准=="是")


CHN_184_s = CHN_184_7_test_results %>% select(...2,uE3检测结果,uE3复测结果,...14,...15) %>% filter(...2 %in% intersect(CHN_184_7_test_results$...2,CHN_184_IE$筛选号))


setdiff(CHN_184_IE$筛选号,CHN_184_s$...2)