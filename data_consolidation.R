install.packages(c("arsenal","compareDF","dplyr","xlsx","readxl","readr","openxlsx"))

library(arsenal)
library(compareDF)
library(dplyr)
library(xlsx)
library(readxl)
library(readr)
library(openxlsx)

# no_value=c("No Value","Cancelled","No result",NA)

##ins align

ins_align <- function(ins_markers,ins_raw,ins_crc){

ins_r = read_csv(ins_raw) %>% 
            filter(`TestName` %in% ins_markers) %>% 
            # filter(! `DoseResult` %in% no_value) %>%
            select(SampleID,DoseResult)

ins_c= read_excel(ins_crc, col_types ="text") %>% 
       rename_at(2,~"DoseResult") %>% 
       mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings() %>%
       rename_at(1,~"SampleID") %>% 
       select(SampleID,DoseResult)

print(setdiff(ins_c$SampleID,ins_r$SampleID))

inter_samples=intersect(ins_r$SampleID,ins_c$SampleID)
compare_df(arrange(ins_r %>% filter(SampleID %in% inter_samples) , SampleID), arrange(ins_c,SampleID)) %>%  create_output_table(output_type = 'xlsx', file_name = paste0(ins_markers[1],"_T_queries.xlsx"))
}


pre_align <- function(pre_markers,pre_raw,pre_crc){

pre_r = read_csv(pre_raw) %>% 
            filter(`Test Name` %in% pre_markers) %>% 
            # filter(! `Result` %in% no_value) %>%
            select(`Sample ID`,Result)

pre_c= read_excel(pre_crc, col_types ="text") %>% 
       rename_at(2,~"Result") %>% 
       mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings() %>%
       rename_at(1,~"Sample ID")%>% 
       select(`Sample ID`,Result)
       

inter_samples=intersect(pre_r$`Sample ID`,pre_c$`Sample ID`)
compare_df(arrange(pre_r %>% filter(`Sample ID` %in% inter_samples) , `Sample ID`), arrange(pre_c, `Sample ID`)) %>%  create_output_table(output_type = 'xlsx', file_name = paste0(pre_markers[1],"_C_queries.xlsx"))
}

ins_label <- function(ins_markers,ins_raw,ins_label){

ins_l= read_excel(ins_label, col_types ="text") %>% 
       rename_at(3,~"Ins_ID") %>% 
       rename_at(1,~"Sample_ID")

ins_r = read_csv(ins_raw, col_types = cols(SampleID = col_character())) %>% 
            filter(`TestName` %in% ins_markers) %>% 
            # filter(! `DoseResult` %in% no_value) %>%
            select(SampleID,DoseResult) %>%
            rename(Ins_ID="SampleID")

print(setdiff(ins_l$Ins_ID,ins_r$Ins_ID))

ins_l %>% left_join(ins_r,by="Ins_ID") %>% select(Sample_ID,DoseResult) %>% arrange(Sample_ID) %>% 
  distinct() %>% 
  group_by(Sample_ID) %>%
  filter(!(n() > 1 & DoseResult=="No result")) %>%
  ungroup() %>%
  write.xlsx(paste0(ins_markers,"_T.xlsx"))
}

pre_label <- function(pre_markers,pre_raw,pre_label){

pre_l= read_excel(pre_label, col_types ="text") %>% 
       rename_at(3,~"Pre_ID") %>% 
       rename_at(1,~"Sample_ID")

pre_r = read_csv(pre_raw, col_types = cols(`Sample ID` = col_character())) %>% 
            filter(`Test Name` %in% pre_markers) %>% 
            # filter(! `Result` %in% no_value) %>%
            select(`Sample ID`,Result) %>%
            rename(Pre_ID="Sample ID")

print(setdiff(pre_l$Pre_ID,pre_r$Pre_ID))

pre_l %>% left_join(pre_r,by="Pre_ID") %>% select(Sample_ID,Result) %>% arrange(Sample_ID) %>% 
  distinct() %>% 
  group_by(Sample_ID) %>%
  filter(!(n() > 1 & Result=="No result")) %>%
  ungroup() %>%
  write.xlsx(paste0(pre_markers,"_C.xlsx"))

}


combine <- function(ss_input,t_input,c_input,out_file,down_l=0,up_l=1000000){

t_mat=read_excel(t_input,col_types ="text") %>% rename_at(1,~"OID") %>% mutate(OID = gsub("R[1-9]","",OID)) %>% 
      mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings() %>%
      select(1,2)

c_mat=read_excel(c_input,col_types ="text")  %>% rename_at(1,~"OID") %>% mutate(OID = gsub("R[1-9]","",OID))%>% 
      mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings() %>%
      select(1,2)

ss_mat=read_excel(ss_input,col_types ="text") %>% rename_at(2,~"OID") 
if(!'备注' %in% names(ss_mat)) ss_mat <- ss_mat %>% mutate(备注 = NA)

print(dim(t_mat))
print(dim(c_mat))
print(dim(ss_mat))


res2 = left_join(ss_mat,t_mat,by="OID")
res3 = left_join(res2,c_mat,by="OID") %>% rename( 唯一可溯源编号 = "OID") %>% 
            mutate(备注 = case_when(
                  ! is.na(备注) ~ 备注,
                  ! .[[ncol(res2)]] %>% as.numeric() %>% between(down_l,up_l) ~ "超出线性范围",
                  ! .[[ncol(res2)+1]] %>% as.numeric() %>% between(down_l,up_l) ~ "超出线性范围",
                  grepl("^>",.[[ncol(res2)]]) ~ "超出线性范围",
                  grepl("^>",.[[ncol(res2)+1]]) ~ "超出线性范围"
                  )) %>% 
            relocate("备注", .after = last_col())

wb <- createWorkbook()
addWorksheet(wb,sheetName = 'Database')
writeData(wb,sheet = 'Database',x = res3)
saveWorkbook(wb,file = paste(out_file,"Database" ,gsub(":", "_", Sys.time( )) ,".xlsx",sep="_"),overwrite = T)

}



combineBT <- function(ss_input,t_input,c_input,b_table,out_file,down_l=0,up_l=1000000){

t_mat=read_excel(t_input,col_types ="text") %>% rename_at(1,~"TID") %>% mutate(TID = gsub("R[1-9]","",TID)) %>% 
      mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings()

c_mat=read_excel(c_input,col_types ="text")  %>% rename_at(1,~"CID") %>% mutate(CID = gsub("R[1-9]","",CID)) %>% 
      mutate_at(2 , ~ifelse(is.na(as.numeric(.)), ., round(as.numeric(.), 3))) %>% suppressWarnings()

ss_mat=read_excel(ss_input, col_types ="text") %>% rename_at(2,~"OID") 
if(!'备注' %in% names(ss_mat)) ss_mat <- ss_mat %>% mutate(备注 = NA)

bt_mat=read_excel(b_table,col_types ="text") %>% rename_at(1,~"OID") %>% rename_at(2,~"TID") %>% rename_at(3,~"CID")

print(dim(t_mat))
print(dim(c_mat))
print(dim(ss_mat))


res2 = ss_mat %>% left_join(bt_mat,by="OID") %>% left_join(t_mat,by="TID") 
                 
res3 = res2 %>% left_join(c_mat,by="CID") %>% 
                  mutate(备注 = case_when(
                  ! is.na(备注) ~ 备注,
                  ! .[[ncol(res2)]] %>% as.numeric() %>% between(down_l,up_l) ~ "超出线性范围",
                  ! .[[ncol(res2)+1]] %>% as.numeric() %>% between(down_l,up_l) ~ "超出线性范围",
                  grepl("^>",.[[ncol(res2)]]) ~ "超出线性范围",
                  grepl("^>",.[[ncol(res2)+1]]) ~ "超出线性范围"
                  )) %>% 
                  select(-c("TID","CID")) %>% rename( 唯一可溯源编号 = "OID")  %>% relocate("备注", .after = last_col())


wb <- createWorkbook()
addWorksheet(wb,sheetName = 'Database')
writeData(wb,sheet = 'Database',x = res3)
saveWorkbook(wb,file = paste(out_file,"Database" ,gsub(":", "_", Sys.time( )) ,".xlsx",sep="_"),overwrite = T)

}

#####merge multiple

library(dplyr)
library(readr)
df <- list.files(path="yourpath", full.names = TRUE) %>% 
  lapply(read_csv(col_types = cols(TestStartDT = col_character(), 
        TestCompleteDT = col_character()))) %>% 
  bind_rows

dim(df)

write_csv(df, "CHN-233_CK MB RI_9000 Results.csv")

####CK-MB

ins_align("CK-MB_IUO","DXI9000 Results_CHN-180_CSV.csv","CHN-180_T_CRA_ZXY_20241230.xlsx")

pre_align("CK-MB","Access2_Results_CHN-180_CSV.csv","CHN-180_C_CRA_ZXY_20241230.xlsx")

combine("CHN-180_SS_CRA_ZXY_20241230.xlsx","CHN-180_T_CRA_ZXY_20241230.xlsx","CHN-180_C_CRA_ZXY_20241230.xlsx","CHN-180")

####RI
ins_align("CK-MB_IUO","DXI9000 Results_CHN-233_CSV.csv","CHN-233_CK-MB_CRA_ZXY_20241230.xlsx" )
ins_align("Myoglobin_IUO","DXI9000 Results_CHN-233_CSV.csv","CHN-233_Myoglobin_CRA_ZXY_20241230.xlsx" )

combine("CHN-233_SS_CRA_ZXY_20241230.xlsx","CHN-233_CK-MB_CRA_ZXY_20241230.xlsx","CHN-233_Myoglobin_CRA_ZXY_20241230.xlsx","CHN-233")


