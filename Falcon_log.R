library(dplyr)
library(tidyverse)
library(xlsx)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(skimr)

parse_log <- function(file){
  log_file_warning = read.csv(file,header  =FALSE, skip = 1,sep = "|") %>% filter(V3 ==" Warning")   
}

####EDC
data_join <- list.files(path = "./warning_logs/", # Identify all CSV files
                       pattern = "*.log", full.names = TRUE) %>% 
  lapply(parse_log) %>%                              # Store all files in list
  reduce(rbind)    