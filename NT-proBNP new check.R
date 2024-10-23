library(dplyr)
library(xlsx)
library(readxl)

###basic check 01

dt = read_excel(dir()[1],sheet = 2,col_names = T)

table(dt$组别)

dt %>% filter(!grepl("^S{1}\\d{6}$", 筛选号))

dtt = dt %>% filter(!grepl("^P{1}\\d{3}$", 样本编号))
table(dtt$样本编号)

table(dt$性别)

sapply(dt$`年龄（岁）`, is.numeric)

dtt = dt %>% filter(!grepl("^\\d+", `NT-proBNP 检测结果（ng/L）`))
table(dtt$`NT-proBNP 检测结果（ng/L）`)

############ 02
dt = read_excel(dir()[1],sheet = 4,col_names = T)

table(dt$组别)

dt %>% filter(!grepl("^S{1}\\d{5}$", 筛选号))

dtt = dt %>% filter(!grepl("^R{1}\\d{5}$", 样本编号))
table(dtt$样本编号)

table(dt$性别)

sapply(dt$`年龄（岁）`, is.numeric)

dtt = dt %>% filter(!grepl("^\\d+", `NT-proBNP 检测结果（ng/L）`))
table(dtt$`NT-proBNP 检测结果（ng/L）`)




#####new

ifNA <- function(x){is.na(x) | x=="NA"}

query <- function(mat){
mat %>% mutate(consec_id = consecutive_id(筛选编号)) %>%
        mutate(dup_No_id = duplicated(筛选编号)) %>%
        mutate(dup_Sample_id = duplicated(样本编号)) %>%
        mutate(consec_id2 = consecutive_id(入组编号)) %>%
        mutate(dup_Enroll_id = duplicated(入组编号)) %>%
        mutate(query = case_when(!grepl("^(A|B|C{1})(\\d{3})",筛选编号) ~ "[筛选编号]的格式不是AXXX、BXXX或CXXX，请核对。",
                                 dup_No_id == TRUE ~ "[筛选编号]有重复值，请核对。",
                                 consec_id!=序号  ~ "[筛选编号]不是连续的数值，请核对。",
                                 入院日期>采样日期 ~ "[入院日期]晚于[采样日期]，请核对。",
                                 # dup_Sample_id == TRUE ~ "[样本编号]有重复值，请核对。",
                                 ! 性别 %in% c("男","女") ~ "填写内容不为<男>或<女>，请核对。",
                                 年龄 < 18 ~ "[年龄]小于<18>岁，请核对。",
                                 # is.na(入院时症状) ~ "数据缺失，请核对。",
                                 # is.na(心超检查时间) ~ "数据缺失，请核对。",
                                 # is.na(心超结果) ~ "数据缺失，请核对。",
                                 # is.na('罗氏NT-proBNP检测时间') | !is.Date('罗氏NT-proBNP检测时间') ~"[罗氏NT-proBNP检测时间]的格式不是YYYY-MM-DD,请核对。"
                                 # is.na('罗氏NT-proBNP检测结果（ng/L）') ~ "数据缺失，请核对。"
                                 # is.na('入院后其它支持诊断的证据（据实选填）') ~ "数据缺失，请核对。",
                                 # is.na('专家确诊诊断') ~ "数据缺失，请核对。",
                                 # !grepl("^(A|B|C{1})(R\\d{3})",入组编号) ~ "[入组编号]的格式不是ARXXX、BRXXX或CRXXX，请核对。",
                                 # dup_Enroll_id == TRUE ~ "[入组编号]有重复值，请核对。",
                                 # consec_id2!=序号  ~ "[入组编号]不是连续的数值，请核对。",
                                 # is.na(入组日期) ~"数据缺失，请核对。",
                                 # is.na('贝克曼NT-proBNP检测时间') | !is.Date('贝克曼NT-proBNP检测时间') ~"[贝克曼NT-proBNP检测时间]的格式不是YYYY-MM-DD,请核对。",
                                 # is.na('贝克曼NT-proBNP检测时间') ~"数据缺失，请核对。",
                                 # 入组日期>'贝克曼NT-proBNP检测时间' ~ "[入组日期]晚于[贝克曼NT-proBNP检测时间]，请核对。",
                                 # is.na('贝克曼NT-proBNP检测结果（ng/L）') ~ "数据缺失，请核对。",
                                 is.na('备注（排除/剔除）原因') ~ "数据缺失，请核对。"
                                 )) %>%
        select (序号,筛选编号,样本编号,入组编号,query)
}

write_excel_csv(read_excel(dir()[1], skip = 2)  %>% query() ,  paste0(dir()[1],"_queries.csv"))

write_excel_csv(read_excel(dir()[2], skip = 2)  %>% query() ,  paste0(dir()[2],"_queries.csv"))
