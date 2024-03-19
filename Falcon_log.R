# grep -s "Warning" C:/Falcon2/warning_logs/{*,.*} >> all_warning.txt


library(dplyr)
library(tidyverse)
library(xlsx)
library(data.table)
library(ggplot2)
library(ggpubr)
library(skimr)



parse_log <- function(file){
  print(file)
  log_file_warning = fread(file,sep = "\t",skip = 1,header = FALSE) %>% 
  extract(V1,c("Time", "Lable","State","Infor"),  "(.+)\\|\\s+(.+)\\|\\s+(.+)\\|\\|\\s+(.+)") %>% 
  filter(State ==" Warning")   
}

####EDC
data_join <- list.files(path = "./warning_logs/", # Identify all CSV files
                       pattern = "*.log", full.names = TRUE) %>% 
  lapply(parse_log) %>%                              # Store all files in list
  reduce(rbind)    