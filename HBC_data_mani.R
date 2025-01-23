library(readxl)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(flextable)
library(janitor)

header_num = 3

data_head <- read_excel(dir()[2],sheet = 2, 
                        n_max = header_num, col_names = FALSE)%>%
             t() %>%       #transpose to a matrix
             as_tibble() %>% fill(V1)  %>% fill(V2) %>% fill(V3) %>%
             mutate(new_names = paste(V1,V2,V3, sep = "_")) %>% 
             mutate(across(new_names, ~str_remove_all(.,"_NA"))) %>% 
             pull(new_names)

dat = readxl::read_excel(dir()[2], sheet = 2,col_names = data_head, skip = header_num) 


total_negative = dat %>% 
        filter(`Routine Test Result_HBsAg...3` == "-") %>%  
        filter(`Routine Test Result_HBsAb...5` == "-") %>% 
        filter(`Routine Test Result_HBeAg...7` == "-") %>% 
        filter(`Routine Test Result_HBeAb...9` == "-") %>% 
        filter(`Routine Test Result_HBV DNA...11` == "NA") %>% 
        select(matches("_Report Result_Interpretation")) 

sapply(total_negative, function(x) prop.table(table(x)))