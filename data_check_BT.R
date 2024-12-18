library(arsenal)
library(compareDF)
library(dplyr)
library(xlsx)
library(readxl)
library(readr)


#######Ferritin

no_value=c("No Value","Cancelled","No result",NA)

bt = read_excel("Gannan_BT.xlsx") 
colnames(bt) = c("Sample_ID","Test_ID","Predicate_ID")

ins_markers=c("Inhibin A_IUO")
ins = read_csv("gannan Dxl9000Results_CSV_100%.csv")%>% 
            filter(! `DoseResult` %in% no_value) %>%
            filter(! `DoseUnit`=="RLU") %>%
            rename(Test_ID="SampleID") %>% 
            mutate(Test_ID = gsub("R[1-9]","",Test_ID))
            
            

predicate_markers=c("InhibinA")
predicate = read_csv("gannan A2Result_CSV_100%.csv",quote = "") %>%  
            filter(! `Result` %in% no_value) %>%
            rename(Predicate_ID = "Sample ID") %>%
#            mutate(Predicate_ID = gsub("R[1-9]","",Predicate_ID))


mat1 = bt %>% 
       left_join(ins  %>% filter(`TestName` %in% ins_markers),by="Test_ID") %>% 
       left_join(predicate  %>% filter(`Test Name` %in% predicate_markers) ,by="Predicate_ID") %>% 
       select(Sample_ID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))


mat2 = read_excel("CHN-216_A_DSF-.xlsx", sheet = 2, skip = 3,col_types ="text") %>% 
       rename_at(7,~"DoseResult") %>% 
       rename_at(8,~"Result") %>% 
       rename(Sample_ID="唯一可溯源编号") %>% 
       select(Sample_ID,DoseResult,Result) %>%
       mutate_at(c('Sample_ID'), as.character)

inter_samples=intersect(mat1$Sample_ID,mat2$Sample_ID)

reform_mat2 =rbind(mat2 %>% distinct(), mat2[duplicated(mat2$Sample_ID),])

create_output_table(compare_df(arrange(mat1 %>% filter(Sample_ID %in% inter_samples) , Sample_ID), reform_mat2 ), output_type = 'xlsx', file_name = paste0(ins_markers[1],"_queries.xlsx"))
