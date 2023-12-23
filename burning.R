library(readxl)
library(mice)
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

select <- dplyr::select

z_scores <- function(data) (ifelse(!is.na(data), abs(data[!is.na(data)]-mean(data[!is.na(data)]))/sd(data[!is.na(data)]),NA ))

model_index <- function(model,marker_name,model_name){
    #index =  model  %>% 
    #model_parameters() %>% 
    #mutate(signif = stars.pval(p))   

    index =  model  %>%
    tidy(conf.int = TRUE) %>% 
    select(c("term","estimate","std.error","conf.low","conf.high","p.value")) %>% 
    mutate(signif = stars.pval(p.value))
    

    main.title <- paste0(marker_name)
    subtitle <- paste0(model_name) %>% strwrap(width = 80) %>% paste(collapse = "\n")
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



## List numerically coded categorical variables
#factorVars <- c("@90天死亡（成活0死亡1）","呕吐（无0有1）")
## Create a variable list. Use dput(names(pbc))
#vars <- c("体重（kg","BMI","TBSA","烧伤指数")

#tableOne <- CreateTableOne(data = empyrosis_t , strata = "呕吐（无0有1）", vars = vars, factorVars = factorVars, smd = TRUE)

test_datasets <- read_excel("test_datasets.xlsx", 
    col_types = c("numeric", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "text", "numeric", "numeric", "numeric", 
        "text", "text", "numeric", "text", 
        "numeric", "text", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "text", "numeric", 
        "text", "text", "text", "text", "text", 
        "numeric", "text", "numeric", "numeric", 
        "text", "text", "numeric", "numeric", 
        "numeric", "text", "text", "text", 
        "numeric", "text", "text", "text"))
#test_datasets %>% select_if(function(x) any(is.na(x))) %>%  %>% mice(method = "pmm")  %>%  complete() %>% as.matrix()
# empyrosis = empyrosis_data2[,6:43]

test_datasets2 <- test_datasets %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(BMI = cut(BMI, breaks = c(-Inf, 28, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(HCTdALB = test_datasets$`@1d_hct`/test_datasets$`@1d_ALB`) %>%
   mutate(`@1d_ALB` = cut(`@1d_ALB`, breaks = c(-Inf, 25, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`24h血糖（mmol/L）` = cut(`24h血糖（mmol/L）`, breaks = c(-Inf, 14, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`= ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==0, "0", "1"))
   



empyrosis_select = test_datasets2 %>% select(-c("性别",matches("@1d")))

dt_f=empyrosis_select %>% select(c("主要结局：耐受0，不耐受1","呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）","年龄","BMI","吸入性损伤...12"))
colnames(dt_f)=c("耐受","呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡","年龄","BMI","吸入性损伤")


row_ha = rowAnnotation(df=as.data.frame(dt_f[,-1]),
                       col=list(呕吐=c('0' = 'blue','1'='red'),
                                呕吐时期=c('0' = 'blue','1'='red','2'='yellow','3'='black'),
                                年龄=c('1'='pink','0'='darkgreen'),
                                BMI=c('1'='pink','0'='darkgreen'),
                                腹泻=c('1'='pink','0'='darkgreen'),
                                CRRT=c('1'='pink','0'='darkgreen'),
                                死亡 =c('1'='pink','0'='darkgreen'),
                                脓毒症=c('1'='pink','0'='darkgreen'),
                                鼻饲=c('1'='pink','0'='darkgreen') ,
                                 吸入性损伤=c('1'='pink','0'='darkgreen')                             
                                # Onset_admission=circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
                                )
                        )

h1 = ComplexHeatmap::Heatmap(empyrosis_select[,c("TBSA","烧伤指数","III度")] %>% mutate_all(as.numeric),
                        column_title = paste0("Key_Value","_Burning"),  
                        right_annotation = row_ha, 
                        left_annotation = rowAnnotation(df=as.data.frame(dt_f[,c('耐受')]),col=list(耐受=c('0' = 'blue','1'='red'))),
                        row_split = dt_f$耐受,
                        col = rev(brewer.pal(10,"RdBu"))
                        )


tableOne1 <- CreateTableOne(vars = colnames(select(test_datasets2, -c("病案号","姓名","体重（kg","主要结局：耐受0，不耐受1","原因","费用","并发症","HCTdALB","吸入性损伤...46"))), 
                           strata = c("主要结局：耐受0，不耐受1"), 
                           data = empyrosis_select)

tb1 = print(tableOne1)     

########index
empyrosis_index = test_datasets2 %>% select(c("主要结局：耐受0，不耐受1",matches("@1d"),"HCTdALB","24h血糖（mmol/L）"))
#colnames(empyrosis_index)=c("耐受","dHB","d_hct","d_plt","d_淋巴细胞","d_TP","d_ALB","d_TBIL","d_DBIL","d_CRE","d_BUN","d_PA","d_LAC","HCTdALB","h血糖")
colnames(empyrosis_index)[1]=c("耐受")

h2= ComplexHeatmap::Heatmap(empyrosis_index %>% select(is.numeric),
                        column_title = paste0("Key_Value","_Blood_index"),  
                        left_annotation = rowAnnotation(df=as.data.frame(empyrosis_index[,c('耐受')]),col=list(耐受=c('0' = 'blue','1'='red'))),
                        row_split = empyrosis_index$耐受,
                        col = rev(brewer.pal(10,"RdBu"))
                        )

tableOne2 <- CreateTableOne(vars = colnames(dplyr::select(empyrosis_index, -c("耐受"))), 
                           strata = c("耐受"), 
                           data = empyrosis_index)

tb2 = print(tableOne2)   



outpdf=paste("Table1","_non.pdf",sep='')
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


t1 <- ggtexttable(as.data.frame(tb1), theme = ttheme("light",base_size = 6,padding = unit(c(4, 4), "mm")))
t2 <- ggtexttable(as.data.frame(tb2), theme = ttheme("light"))      

print(t1)
print(t2)

dev.off()

##regression model

dt_select = test_datasets2 %>%  
select(-c("病案号","姓名","体重（kg","身高m","天数","中性粒细胞","吸入性损伤...46"))   %>% rename_with(~str_remove(., '[@）]+')) %>% mutate_if(is.character, as.factor) %>% na.omit()
colnames(dt_select) <- gsub('[(（].*','',colnames(dt_select))
colnames(dt_select) [25] = "耐受"
#colnames(dt_select) [40] = "吸入性损伤"
colnames(dt_select) [8] = "吸入性损伤"

biomarker="ALL_X"

#dt_select = empyrosis_data2 %>% 
#            select(c("主要结局：耐受0，不耐受1","年龄","性别","BMI","TBSA","烧伤指数","III度",matches("@1d"),"呕吐（无0有1）","呕吐（无呕吐0，休克期1，非休克期2，两者都有3）","脓毒症（无0有1）","鼻饲（有1无0）","腹泻","CRRT","@90天死亡（成活0死亡1）")) %>% 
#            mutate_at(vars(!c("BMI","TBSA","烧伤指数","III度",matches("@1d"))), as.factor) %>%
#            mutate_at(vars("年龄"), as.numeric)


# index <- sample(nrow(dt_na_removed),nrow(dt_na_removed)*0.80)
# dt_train = dt_na_removed[index,]
# dt_test = dt_na_removed[-index,]



fit.glm = glm(耐受 ~ 年龄+性别+BMI+吸入性损伤+天数+基础疾病+`24h血糖`+TBSA+烧伤指数+III度+脓毒症+HCTdALB+`1d_ALB`+`1dHB`+`1d_plt`+`1d_淋巴细胞`+`1d_TP`+`1d_TBIL`+`1d_DBIL`+`1d_CRE`+`1d_BUN`+`1d_PA`+`1d_LAC`, family = binomial,data=dt_select) 

fit.glm.select =  fit.glm %>% select_parameters()

t.glm <- model_index(fit.glm,biomarker,"Generalized Linear Model")

t.glm.select <- model_index(fit.glm.select ,biomarker,"Generalized Linear Model with Select variable")


######

biomarker="Select_X & remove random effects"

fit.mem <- glmer(耐受 ~ (1 | 原因) + (1 | 启动时间) +年龄+性别+BMI+吸入性损伤+天数+基础疾病+`24h血糖`+TBSA+烧伤指数+III度+脓毒症+HCTdALB+`1d_ALB`+`1dHB`+`1d_plt`+`1d_淋巴细胞`+`1d_TP`+`1d_TBIL`+`1d_DBIL`+`1d_CRE`+`1d_BUN`+`1d_PA`+`1d_LAC`, data = dt_select, family = binomial, control = glmerControl(optimizer = "bobyqa")) 

fit.mem.select <- glmer(耐受 ~ (1 | 原因) + (1 | 启动时间)+ 年龄 + 吸入性损伤 + 烧伤指数 + HCTdALB+ `1d_plt`+`1d_TBIL`, data = dt_select, family = binomial)

t.mem <- model_index(fit.mem,biomarker,"Mix Model random effects removed")

t.mem.select <- model_index(fit.mem.select ,biomarker,"Mix Model with Select variable")

# ROC curve for first model

validate_datasets <- read_excel("validate_datasets.xlsx")

dt_test= validate_datasets %>%
   mutate(HCTdALB = validate_datasets$`@1d_hct`/validate_datasets$`@1d_ALB`) %>% 
   select(c("主要结局：耐受0，不耐受1","原因","启动时间（天）","年龄","吸入性损伤","烧伤指数","HCTdALB","@1d_plt","@1d_TBIL")) %>%
   mutate(年龄 = str_remove(年龄, "岁")) %>% mutate_at(vars("年龄"), as.numeric) %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>% 
   mutate_at(vars(c("主要结局：耐受0，不耐受1","原因","启动时间（天）","年龄","吸入性损伤")), as.factor) %>% 
   rename_with(~str_remove(., '[@）]+'))  %>% na.omit()
colnames(dt_test)[1:3] = c("耐受", "原因","启动时间")


glmmod_plot <- ROC_curve(fit.glm,dt_select,dt_select$耐受,'广义线性模型_未筛选因子_训练集') 
glm.select.mod_plot <- ROC_curve(fit.glm.select,dt_select,dt_select$耐受,'广义线性模型_筛选因子_训练集') 
glm.select.mod_plot_test <- ROC_curve(fit.glm.select,dt_test,dt_test$耐受,'广义线性模型_筛选因子_测试集',test='T') 

# ROC curve for second model

memmod_plot <- ROC_curve(fit.mem,dt_select,dt_select$耐受,'混合模型去除随机效应_未筛选因子_训练集') 
mem.select.mod_plot <- ROC_curve(fit.mem.select ,dt_select,dt_select$耐受,'混合模型去除随机效应_筛选因子_训练集') 
mem.select.mod_plot_test <- ROC_curve(fit.mem.select ,dt_test,dt_test$耐受,'混合模型去除随机效应_筛选因子_测试集',test='T') 

outpdf=paste("resgression","_profile_new.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

pcom.glm = ggarrange(t.glm,t.glm.select , ncol = 2, nrow = 1,widths=c(1, 1))
print(pcom.glm)
ggpubr::ggarrange(glmmod_plot, glm.select.mod_plot, nrow = 1)


pcom.mem = ggarrange(t.mem,t.mem.select , ncol = 2, nrow = 1,widths=c(1, 1))
print(pcom.mem)
ggpubr::ggarrange(memmod_plot, mem.select.mod_plot, nrow = 1)


ggpubr::ggarrange(glm.select.mod_plot_test, mem.select.mod_plot_test, nrow = 1)


dev.off()


##################
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3fselect)
library(mlr3viz)
library(skimr)

dt_train = dt_select  %>% 
    select(c("耐受","年龄","性别","BMI","吸入性损伤","天数","基础疾病","24h血糖","TBSA","烧伤指数","III度","脓毒症","HCTdALB","1d_ALB","1dHB","1d_plt","1d_淋巴细胞","1d_TP","1d_TBIL","1d_DBIL","1d_CRE","1d_BUN","1d_PA","1d_LAC"))

colnames(dt_train)=c("Res","Age","Sex","Bmi","Inhalation_injury","Day","Underlying_disease","Glu_24h","TBSA","Burn_index","III_index","Sepsis","HCTdALB","ALB","HB","Plt","Lymphocyte","TP","TBIL","DBIL","CRE","BUN","PA","LAC")

dt_test = validate_datasets %>%
   mutate(年龄 = str_remove(年龄, "岁")) %>% mutate_at(vars("年龄"), as.numeric) %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(BMI = cut(BMI, breaks = c(-Inf, 28, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(HCTdALB = validate_datasets$`@1d_hct`/validate_datasets$`@1d_ALB`) %>% 
   mutate(`@1d_ALB` = cut(`@1d_ALB`, breaks = c(-Inf, 25, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`24h血糖` = cut(`24h血糖`, breaks = c(-Inf, 14, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`= ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==0, "0", "1")) %>%
   mutate(天数 = cut(天数, breaks = c(-Inf, 4, Inf), right = FALSE, labels = c("0", "1"))) %>%
   select(c("主要结局：耐受0，不耐受1","年龄","性别","BMI","吸入性损伤","天数","基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)","24h血糖","TBSA","烧伤指数","III度","脓毒症（无0有1）","HCTdALB","@1d_ALB","@1dHB","@1d_plt","@1d_淋巴细胞","@1d_TP","@1d_TBIL","@1d_DBIL","@1d_CRE","@1d_BUN","@1d_PA","@1d_LAC")) %>% 
   mutate_if(is.character, as.factor) %>% na.omit()
colnames(dt_test)=c("Res","Age","Sex","Bmi","Inhalation_injury","Day","Underlying_disease","Glu_24h","TBSA","Burn_index","III_index","Sepsis","HCTdALB","ALB","HB","Plt","Lymphocyte","TP","TBIL","DBIL","CRE","BUN","PA","LAC")

skimr::skim(dt_test)

dataset = rbind(dt_train,dt_test) %>%
   select(c("Res", "Age","Inhalation_injury","Burn_index","HCTdALB","Plt","TBIL"))


task = TaskClassif$new("shao_test", dataset, target = "Res")

custom = rsmp("custom")
train_sets = list(1:nrow(dt_train), (nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)))
test_sets = list((nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)), 1:nrow(dt_train))
custom$instantiate(task, train_sets, test_sets)

custom1 = rsmp("custom")
train_sets = list(1:nrow(dt_train))
test_sets = list((nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)))
custom1$instantiate(task, train_sets, test_sets)

#task = TaskClassif$new("shao_test", dt_train, target = "Res")


at_log_reg = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.log_reg", predict_type = "prob"),resampling = rsmp("cv", folds = 3), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_ranger = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.ranger", predict_type = "prob"),resampling = rsmp("cv", folds = 3), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
# at_svm = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.svm", predict_type = "prob"),resampling = rsmp("cv", folds = 5), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_rpart = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.rpart", predict_type = "prob"),resampling = rsmp("cv", folds = 3), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)

auto <- auto_fselector(
  fselector = fs("random_search"),
  learner = lrn("classif.log_reg", predict_type = "prob"),
  resampling = rsmp("holdout"),
  measure = msr("classif.auc"),
  terminator = trm("evals", n_evals = 10)
)

learners <- c(at_log_reg,  at_ranger,at_rpart,auto)
measures <- msrs(c("classif.auc", "classif.bacc", "classif.bbrier"))

#Benchmarking
set.seed(372)
design = benchmark_grid(tasks =task, learners = learners, resamplings = rsmp("cv", folds = 5) )
bmr = benchmark(design, store_models = TRUE)
results <- bmr$aggregate(measures)
print(results)
autoplot(bmr, measure = msr("classif.auc"))
autoplot(bmr, type = "roc")

#######################################

design = benchmark_grid(
  tasks = task,
  learners = lrns(c("classif.log_reg", "classif.ranger","classif.rpart"),predict_type = "prob"),
  # resampling = rsmp("cv", folds = 5)
  resamplings = custom 
)


bmr <-  benchmark(design)

autoplot(bmr, type = "roc")

bmr$score(msr("classif.auc"))

aggr <-  bmr$aggregate(msrs(c("classif.acc", "time_train")))
as.data.table(aggr)[, .(learner_id, classif.acc, time_train)]











####################backup code









learner_logreg = lrn("classif.log_reg")
learner_logreg$train(task)

######80 training sets
train_set = sample(task$row_ids, 0.8 * task$nrow)
test_set = setdiff(task$row_ids, train_set)
learner_logreg$train(task, row_ids = train_set)

summary(learner_logreg$model) 


auto <- auto_fselector(
  fselector = fs("random_search"),
  learner = lrn("classif.log_reg"),
  resampling = rsmp("holdout"),
  measure = msr("classif.acc"),
  terminator = trm("evals", n_evals = 10)
)
grid <-  benchmark_grid(
  task = task,
  learner = list(auto, lrn("classif.log_reg")),
  resampling = rsmp("cv", folds = 3)
)

bmr <-  benchmark(grid)

aggr <-  bmr$aggregate(msrs(c("classif.acc", "time_train")))
as.data.table(aggr)[, .(learner_id, classif.acc, time_train)]


#####random forest
learner_rf = lrn("classif.ranger", importance = "permutation")
learner_rf$train(task, row_ids = train_set)

importance = as.data.table(learner_rf$importance(), keep.rownames = TRUE)
colnames(importance) = c("Feature", "Importance")
ggplot(data=importance,aes(x = reorder(Feature, Importance), y = Importance)) + geom_col() + coord_flip() + xlab("")


lrn_ranger = lrn("classif.ranger", predict_type = "prob")
splits = mlr3::partition(task, ratio = 0.8)

lrn_ranger$train(task, splits$train)
prediction = lrn_ranger$predict(task, splits$test)
autoplot(prediction, type = "roc")
prediction$score(msr("classif.auc"))



library(patchwork)

design = benchmark_grid(
  tasks = task,
  learners = lrns(c("classif.log_reg", "classif.ranger"),predict_type = "prob"),
  resamplings = rsmp("cv", folds = 5)
)
bmr = benchmark(design)
autoplot(bmr, type = "roc") +
  plot_layout(guides = "collect")


autoplot(bmr, measure = msr("classif.auc"))
autoplot(bmr, type = "roc")


rr$learners[[1]]$predict_newdata(dat_1[rr$resampling$test_set(1), ], task = rr$task)