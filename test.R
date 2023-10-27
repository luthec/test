library(readxl)
library(tableone)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridGraphics)
library(tidyverse)
library(lmerTest)
library(lme4)
library(parameters)
library(broom)
library(broom.mixed)
library(gtools)
library(pROC)
library(ggplot2)
library(report)
#library(cvms)
#library(performance)
#library(Cairo)

rm(list = ls())

model_index <- function(model,marker_name,model_name){
    index =  model  %>% 
    model_parameters() %>% 
    mutate(signif = stars.pval(p))   

    main.title <- paste0(marker_name)
    subtitle <- paste0(model_name) %>%
    strwrap(width = 80) %>%
    paste(collapse = "\n")
    ggtexttable(index, theme = ttheme("light")) %>%
    tab_add_title(text = subtitle, face = "plain", size = 10) %>%
    tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) %>%
    tab_add_footnote(text = "* means significant levels", size = 10, face = "italic")

}


ROC_curve <- function(model,data_mat,y,model_name,test=NA){
  if (is.na(test)){
    probs_glmmod <- predict(model, type = 'response')
    roc_glmmod  <- roc(response = y, predictor = probs_glmmod)
  }else {
    probs_glmmod <- predict(model, newdata = data_mat, type="response",allow.new.levels = TRUE)
    roc_glmmod  <- roc(response = y, predictor = probs_glmmod)
  }

    ggroc(roc_glmmod, legacy.axes = TRUE) +
    labs(x = '假阳性率', y = '真阳性率',
    title = paste0('ROC curve\n Y~',model_name)) +
    annotate('text', x = .5, y = .5, label = paste0('AUC: ', round(auc(roc_glmmod), digits = 2)))
}


empyrosis_data2 <- read_excel("empyrosis_data2.xlsx", 
    col_types = c("text", "text", "text", 
         "numeric", "numeric", "numeric", 
         "text", "numeric", "numeric", "numeric", 
         "text", "numeric", "numeric", "numeric", 
         "numeric", "text", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "text", "text", 
         "text", "text", "text", "text", "text", 
         "numeric", "text", "numeric", "numeric", 
        "text", "text", "numeric", "numeric", 
         "text"))

## List numerically coded categorical variables
#factorVars <- c("@90天死亡（成活0死亡1）","呕吐（无0有1）")
## Create a variable list. Use dput(names(pbc))
#vars <- c("体重（kg","BMI","TBSA","烧伤指数")

#tableOne <- CreateTableOne(data = empyrosis_t , strata = "呕吐（无0有1）", vars = vars, factorVars = factorVars, smd = TRUE)

empyrosis = empyrosis_data2[,6:43]

empyrosis_select = empyrosis %>% select(-c("性别",matches("@1d")))

dt_f=empyrosis_select %>% select(c("主要结局：耐受0，不耐受1","呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）"))
colnames(dt_f)=c("耐受","呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡")

row_ha = rowAnnotation(df=as.data.frame(dt_f[,-1]),
                       col=list(呕吐=c('0' = 'blue','1'='red'),
                                呕吐时期=c('0' = 'blue','1'='red','2'='yellow','3'='black'),
                                脓毒症=c('1'='pink','0'='darkgreen'),
                                鼻饲=c('1'='pink','0'='darkgreen'),
                                腹泻=c('1'='pink','0'='darkgreen'),
                                CRRT=c('1'='pink','0'='darkgreen'),
                                死亡 =c('1'='pink','0'='darkgreen')                              
                                # Onset_admission=circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
                                )
                        )

h1 = ComplexHeatmap::Heatmap(empyrosis_select[,3:5],
                        column_title = paste0("Key_Value","_Burning"),  
                        right_annotation = row_ha, 
                        left_annotation = rowAnnotation(df=as.data.frame(dt_f[,c('耐受')]),col=list(耐受=c('0' = 'blue','1'='red'))),
                        row_split = dt_f$耐受,
                        col = rev(brewer.pal(10,"RdBu"))
                        )


tableOne <- CreateTableOne(vars = colnames(select(empyrosis_select, -c("主要结局：耐受0，不耐受1","原因","费用","并发症"))), 
                           strata = c("主要结局：耐受0，不耐受1"), 
                           data = empyrosis_select)

tb1 = print(
  tableOne,
  nonnormal = c("TBSA","烧伤指数"),exact = c("@90天死亡（成活0死亡1）"),
  showAllLevels = TRUE)     

########index
empyrosis_index = empyrosis %>% select(c("主要结局：耐受0，不耐受1",matches("@1d")))
colnames(empyrosis_index)[1]="耐受"

h2= ComplexHeatmap::Heatmap(empyrosis_index[,-1],
                        column_title = paste0("Key_Value","_Blood_index"),  
                        left_annotation = rowAnnotation(df=as.data.frame(empyrosis_index[,c('耐受')]),col=list(耐受=c('0' = 'blue','1'='red'))),
                        row_split = empyrosis_index$耐受,
                        col = rev(brewer.pal(10,"RdBu"))
                        )

tableOne <- CreateTableOne(vars = colnames(select(empyrosis_index, -c("耐受"))), 
                           strata = c("耐受"), 
                           data = empyrosis_index)

tb2 = print(
  tableOne,
  nonnormal = c("TBSA","烧伤指数"),
  showAllLevels = TRUE)   

##regression model

biomarker="ALL_X"

dt_select = empyrosis_data2 %>% 
            select(c("主要结局：耐受0，不耐受1","年龄","性别","BMI","TBSA","烧伤指数","III度",matches("@1d"),"呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）")) %>% 
            mutate_at(vars(!c("BMI","TBSA","烧伤指数","III度",matches("@1d"))), as.factor) %>%
            mutate_at(vars("年龄"), as.numeric)


colnames(dt_select)[c(1,20:26)]=c("耐受","呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡")

dt_na_removed <- na.omit(dt_select)

index <- sample(nrow(dt_na_removed),nrow(dt_na_removed)*0.80)
dt_train = dt_na_removed[index,]
dt_test = dt_na_removed[-index,]

fit.glm = glm(耐受 ~ 年龄+性别+BMI+TBSA+烧伤指数+III度+脓毒症+`@1d_TBIL`+`@1d_DBIL`+`@1d_BUN`+`@1d_LAC`+`@1d_CRE`+`@1d_hct`+`@1d_ALB`+`@1dHB`+`@1d_淋巴细胞`+`@1d_plt`+`@1d_PA`+`@1d_TP`, family = binomial,data=dt_train) 

fit.glm.select =  fit.glm %>% select_parameters()

t.glm <- model_index(fit.glm,biomarker,"Generalized Linear Model")

t.glm.select <- model_index(fit.glm.select ,biomarker,"Generalized Linear Model with Select")

# ROC curve for first model

glmmod_plot <- ROC_curve(fit.glm,dt_train,dt_train$耐受,'未筛选因子') 
glmmod_plot_test <- ROC_curve(fit.glm,dt_test,dt_test$耐受,'未筛选因子',test='T') 

######

biomarker="Select_X & remove random effects"

dt_select2 = empyrosis_data2 %>% 
            select(c("原因","受伤至入院时间（h）","启动时间（天）","主要结局：耐受0，不耐受1","年龄","性别","BMI","TBSA","烧伤指数","III度",matches("@1d"),"呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）")) %>% 
            mutate_at(vars(!c("BMI","TBSA","烧伤指数","III度",matches("@1d"))), as.factor) %>%
            mutate_at(vars("年龄"), as.numeric) %>%
            mutate_at(vars(c("原因","受伤至入院时间（h）","启动时间（天）")), as.factor)


colnames(dt_select2)[c(1:4,23:29)]=c("原因","入院时间","启动时间","耐受","呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡")

dt_na_removed <- na.omit(dt_select2)


index <- sample(nrow(dt_na_removed),nrow(dt_na_removed)*0.80)
dt_train = dt_na_removed[index,]
dt_test = dt_na_removed[-index,]


fit.mem <- glmer(耐受 ~ (1 | 原因) +(1 | 入院时间) +(1 | 启动时间) +年龄+性别+BMI+TBSA+烧伤指数+III度+脓毒症+`@1d_TBIL`+`@1d_DBIL`+`@1d_BUN`+`@1d_LAC`+`@1d_CRE`+`@1d_hct`+`@1d_ALB`+`@1dHB`+`@1d_淋巴细胞`+`@1d_plt`+`@1d_PA`+`@1d_TP`, data = dt_train, family = binomial, control = glmerControl(optimizer = "bobyqa")) 
fit.mem2 <- glmer(耐受 ~ (1 | 原因) +(1 | 入院时间) +(1 | 启动时间) +年龄+性别+BMI+TBSA+烧伤指数+III度+脓毒症+`@1d_TBIL`, data = dt_train, family = binomial) %>%
  select_parameters()



t.mem <- model_index(fit.mem,biomarker,"Generalized Linear Model")

# ROC curve for second model

memmod_plot <- ROC_curve(fit.mem,dt_train,dt_train$耐受,'去除随机效应') 
memmod_plot_test <- ROC_curve(fit.mem,dt_test,dt_test$耐受,'去除随机效应',test='T') 



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
t2 <- ggtexttable(as.data.frame(tb2), theme = ttheme("light"))      
tcom = ggarrange(t1,t2, ncol = 2, nrow = 1,widths=c(1, 1))
print(tcom)

pcom = ggarrange(t.glm,t.mem, ncol = 2, nrow = 1,widths=c(1, 1))
print(pcom)

ggpubr::ggarrange(glmmod_plot, memmod_plot, nrow = 1)

ggpubr::ggarrange(glmmod_plot_test, memmod_plot_test, nrow = 1)

dev.off()
