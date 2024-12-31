library(arsenal)
library(compareDF)
library(dplyr)
library(xlsx)
library(readxl)
library(readr)
library(openxlsx)

no_value=c("No Value","Cancelled","No result",NA)

##ins align

ins_align <- function(ins_markers,ins_raw,ins_crc){

ins_r = read_csv(ins_raw) %>% 
            filter(`TestName` %in% ins_markers) %>% 
            filter(! `DoseResult` %in% no_value) %>%
            select(SampleID,DoseResult)

ins_c= read_excel(ins_crc, col_types ="text") %>% 
       rename_at(2,~"DoseResult") %>% 
       rename_at(1,~"SampleID")

inter_samples=intersect(ins_r$SampleID,ins_c$SampleID)
compare_df(arrange(ins_r %>% filter(SampleID %in% inter_samples) , SampleID), ins_c) %>%  create_output_table(output_type = 'xlsx', file_name = paste0(ins_markers[1],"_T_queries.xlsx"))
}


pre_align <- function(pre_markers,pre_raw,pre_crc){

pre_r = read_csv(pre_raw) %>% 
            filter(`Test Name` %in% pre_markers) %>% 
            filter(! `Result` %in% no_value) %>%
            select(`Sample ID`,Result)

pre_c= read_excel(pre_crc, col_types ="text") %>% 
       rename_at(2,~"Result") %>% 
       rename_at(1,~"Sample ID")

inter_samples=intersect(pre_r$`Sample ID`,pre_c$`Sample ID`)
compare_df(arrange(pre_r %>% filter(`Sample ID` %in% inter_samples) , `Sample ID`), pre_c) %>%  create_output_table(output_type = 'xlsx', file_name = paste0(pre_markers[1],"_C_queries.xlsx"))
}



combine <- function(ss_input,t_input,c_input,out_file){

t_mat=read_excel(t_input,sheet = 1) %>% rename_at(1,~"OID")

c_mat=read_excel(c_input,sheet = 1)  %>% rename_at(1,~"OID") %>% mutate(OID = gsub("R[1-9]","",OID))

ss_mat=read_excel(ss_input,sheet = 1) %>% rename_at(2,~"OID") %>% mutate(OID = gsub("R[1-9]","",OID))

print(dim(t_mat))
print(dim(c_mat))
print(dim(ss_mat))


res2 = left_join(ss_mat,t_mat,by="OID")
res3 = left_join(res2,c_mat,by="OID") %>% rename( 唯一可溯源编号 = "OID")


print(dim(res2))
print(dim(res3))

wb <- createWorkbook()
addWorksheet(wb,sheetName = 'Database')
writeData(wb,sheet = 'Database',x = res3)
saveWorkbook(wb,file = paste(out_file,"Database" ,gsub(":", "_", Sys.time( )) ,".xlsx",sep="_"),overwrite = T)

}



####CK-MB


ins_align("CK-MB_IUO","DXI9000 Results_CHN-180_CSV.csv","CHN-180_T_CRA_ZXY_20241230.xlsx")

pre_align("CK-MB","Access2_Results_CHN-180_CSV.csv","CHN-180_C_CRA_ZXY_20241230.xlsx")


