library(readxl)
library(tableone)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridGraphics)
#library(Cairo)

empyrosis_data1 <- read_excel("empyrosis_data1.xlsx", 
     col_types = c("text", "text", "numeric", 
         "numeric", "numeric", "text", "numeric", 
         "numeric", "numeric", "text", "numeric", 
         "numeric", "text", "numeric", "text", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "text", "text", "text", 
         "text", "numeric", "text", "numeric", 
         "numeric", "text", "text", "numeric", 
         "numeric"))


empyrosis_t =  empyrosis_data1 %>% select(c("体重（kg","BMI","TBSA","烧伤指数","@90天死亡（成活0死亡1）","呕吐（无0有1）"))

## List numerically coded categorical variables
#factorVars <- c("@90天死亡（成活0死亡1）","呕吐（无0有1）")
## Create a variable list. Use dput(names(pbc))
#vars <- c("体重（kg","BMI","TBSA","烧伤指数")

#tableOne <- CreateTableOne(data = empyrosis_t , strata = "呕吐（无0有1）", vars = vars, factorVars = factorVars, smd = TRUE)


empyrosis_select = empyrosis_data1[,7:40] %>% select(-c("性别",matches("@1d")))

dt_f=empyrosis_select %>% select(c("呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）"))
colnames(dt_f)=c("呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡")

row_ha = rowAnnotation(df=as.data.frame(dt_f[,-1]),
                       col=list(呕吐时期=c('0' = 'blue','1'='red','2'='yellow','3'='black'),
                                脓毒症=c('1'='pink','0'='darkgreen'),
                                鼻饲=c('1'='pink','0'='darkgreen'),
                                腹泻=c('1'='pink','0'='darkgreen'),
                                CRRT=c('1'='pink','0'='darkgreen'),
                                死亡 =c('1'='pink','0'='darkgreen')                              
                                # Onset_admission=circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
                                )
                        )

h1 = ComplexHeatmap::Heatmap(na.omit(empyrosis_select[,1:3]),
                        column_title = paste0("Key_Value","_Burning"),  
                        right_annotation = row_ha, 
                        left_annotation = rowAnnotation(df=as.data.frame(dt_f[,c('呕吐')]),col=list(呕吐=c('0' = 'blue','1'='red'))),
                        row_split = dt_f$呕吐,
                        col = rev(brewer.pal(10,"RdBu"))
                        )


tableOne <- CreateTableOne(vars = colnames(select(empyrosis_select, -c("呕吐（无0有1）","费用","并发症"))), 
                           strata = c("呕吐（无0有1）"), 
                           data = empyrosis_select)

tb1 = print(
  tableOne,
  nonnormal = c("TBSA","烧伤指数"),exact = c("@90天死亡（成活0死亡1）"),
  showAllLevels = TRUE)     

########index
empyrosis_index = empyrosis_data1[,7:40] %>% select(c("呕吐（无0有1）",matches("@1d")))
colnames(empyrosis_index)[1]="呕吐"

h2= ComplexHeatmap::Heatmap(empyrosis_index[,-1],
                        column_title = paste0("Key_Value","_Blood_index"),  
                        left_annotation = rowAnnotation(df=as.data.frame(empyrosis_index[,c('呕吐')]),col=list(呕吐=c('0' = 'blue','1'='red'))),
                        row_split = empyrosis_index$呕吐,
                        col = rev(brewer.pal(10,"RdBu"))
                        )

tableOne <- CreateTableOne(vars = colnames(select(empyrosis_index, -c("呕吐"))), 
                           strata = c("呕吐"), 
                           data = empyrosis_index)

tb2 = print(
  tableOne,
  nonnormal = c("TBSA","烧伤指数"),
  showAllLevels = TRUE)   

outpdf=paste("res","_profile_new.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h1,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h2,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)



t1 <- ggtexttable(as.data.frame(tb1), theme = ttheme("light"))      
tcom = ggarrange(t1,t2, ncol = 2, nrow = 1,widths=c(1, 1))
print(pcom)

dev.off()
