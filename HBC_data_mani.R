library(readxl)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(flextable)
library(janitor)

header_num = 3

data_head <- read_excel(dir()[1],sheet = 2, 
                        n_max = header_num, col_names = FALSE)%>%
             t() %>%       #transpose to a matrix
             as_tibble() %>% fill(V1)  %>% fill(V2) %>% fill(V3) %>%
             mutate(new_names = paste(V1,V2,V3, sep = "_")) %>% 
             mutate(across(new_names, ~str_remove_all(.,"_NA"))) %>% 
             pull(new_names)

dat = readxl::read_excel(dir()[1], sheet = 2,col_names = data_head, skip = header_num) 

########select 2 validation
dat  %>% filter(`BEC DxI 9000_Report Result_Interpretation` == "Reactive") %>%
         filter(`Abbott Alinity_Report Result_Interpretation` == "Nonreactive") %>%
         select(matches("_Report Result_Interpretation")) %>%
         rename_all(function(x) gsub("_.+$", "", x)) %>%
         select("BEC DxI 9000","Abbott Alinity","Siemens Atellica")


dat  %>% filter(`BEC DxI 9000_Report Result_Interpretation` == "Nonreactive") %>%
         filter(`Abbott Alinity_Report Result_Interpretation` == "Reactive") %>%
         select(matches("_Report Result_Interpretation")) %>%
         rename_all(function(x) gsub("_.+$", "", x)) %>%
         select("BEC DxI 9000","Abbott Alinity","Siemens Atellica")



####total_negative
total_negative = dat %>% 
        filter(`Routine Test Result_HBsAg...3` == "-") %>%  
        filter(`Routine Test Result_HBsAb...5` == "-") %>% 
        filter(`Routine Test Result_HBeAg...7` == "-") %>% 
        filter(`Routine Test Result_HBeAb...9` == "-") %>% 
        filter(`Routine Test Result_HBV DNA...11` == "NA") %>% 
        select(matches("_Report Result_Interpretation"))  %>%
        rename_all(function(x) gsub("_.+$", "", x))

main.title="Never infected and no evidence of immunization"
tn <- sapply(total_negative, function(x) table(x)) %>% t() %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))
tn0 <- sapply(total_negative, function(x) prop.table(table(x))) %>% t() %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))
tn_f <- total_negative %>% filter(`BEC DxI 9000` == "Reactive") %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))


vaccination_negative = dat %>% 
        filter(`Routine Test Result_HBsAg...3` == "-") %>%  
        filter(`Routine Test Result_HBsAb...5` == "+") %>% 
        filter(`Routine Test Result_HBeAg...7` == "-") %>% 
        filter(`Routine Test Result_HBeAb...9` == "-") %>% 
        filter(`Routine Test Result_HBV DNA...11` == "NA") %>% 
        select(matches("_Report Result_Interpretation"))  %>%
        rename_all(function(x) gsub("_.+$", "", x))

main.title="Immune due to hepatitis B vaccination"
vn <- sapply(vaccination_negative, function(x) table(x)) %>% t() %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))
vn0 <- sapply(vaccination_negative, function(x) prop.table(table(x))) %>% t() %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))
vn_f <- vaccination_negative %>% filter(`BEC DxI 9000` == "Reactive") %>% ggtexttable(theme = ttheme("light")) %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))

outpdf=paste("Table2","presume_negative_samples.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

print(tn)
print(tn0)
print(tn_f)

print(vn)
print(vn0)
print(vn_f)

dev.off()