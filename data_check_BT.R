#######Ferritin

ins_markers=c("Inhibin A_IUO")
ins = read_csv(dir()[2])%>% 
            rename(Test_ID="SampleID")

predicate_markers=c("InhibinA")
predicate = read_csv(dir()[1],quote = "") %>% 
            rename(Predicate_ID="Sample ID")

bt = read_excel(dir()[3]) 
colnames(bt) = c("Sample_ID","Test_ID","Predicate_ID")

mat1 = bt %>% 
       left_join(ins  %>% filter(`TestName` %in% ins_markers),by="Test_ID") %>% 
       left_join(predicate  %>% filter(`Test Name` %in% predicate_markers) 
                            %>% filter(`Test Name` %in% predicate_markers)
       ,by="Predicate_ID") %>% 
       select(Sample_ID,DoseResult,Result)

# write_excel_csv(ins %>% filter(TestName %in% markers) %>% select(SampleID,TestName,DoseResult),  paste0(markers[1],"_wrongnumbers.csv"))


mat2 = read_excel(dir()[4], sheet = 1, skip = 3,col_types ="text") %>% 
       rename_at(7,~"DoseResult") %>% 
       rename_at(8,~"Result") %>% 
       rename(Sample_ID="唯一可溯源编号") %>% 
       select(Sample_ID,DoseResult,Result) %>%
       mutate_at(c('Sample_ID'), as.character)

inter_samples=intersect(mat1$Sample_ID,mat2$Sample_ID)

create_output_table(compare_df(arrange(mat1 %>% filter(Sample_ID %in% inter_samples) , Sample_ID),arrange(mat2 %>% filter(Sample_ID %in% inter_samples) , Sample_ID)), output_type = 'xlsx', file_name = paste0(markers[1],"_queries.xlsx"))
