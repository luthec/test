# grep -s "Warning" C:/Falcon2/warning_logs/{*,.*} >> all_warning.txt
# grep -s "Error" all_warning.txt > error_log.txt

library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(DataExplorer)
library(readr)



parse_log <- function(file){
  print(file)
  log_file_warning = fread(file,sep = "\t",skip = 1,header = FALSE) %>% 
  extract(V1,c("Filename", "Time", "Lable","State","Infor"),  "(.+)\\.log\\:(.+)\\|\\s+(.+)\\|\\s+(.+)\\|\\|\\s+(.+)") %>% 
  filter(State ==" Warning")   
}

####EDC
data_join <- list.files(path = "./warning_logs/", # Identify all CSV files
                       pattern = "*.log", full.names = TRUE) %>% 
  lapply(parse_log) %>%                              # Store all files in list
  reduce(rbind)    


log_file_error = fread("error_log.txt",sep = "\t",skip = 1,header = FALSE) %>% 
    extract(V1,c("Filename", "Time", "Lable","State","Infor"),  "(.+)\\.log\\:(.+)\\|\\s+(.+)\\|\\s+(.+)\\|\\|\\s+(.+)")


log_file_error2 = log_file_error %>% 
  extract(Infor, c("Trace","Class","Error_Infro","Number","Error","Detail","Time2"),"(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)") 

create_report(log_file_error2)

