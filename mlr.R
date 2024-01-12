##################

library(readxl)
library(mice)
library(tableone)
library(dplyr)
library(ggpubr)

library(gridGraphics)
library(tidyverse)

library(pROC)
library(ggplot2)
library(mlr3)
library(mlr3verse)
library(mlr3learners)
library(mlr3tuning)
library(mlr3fselect)
library(mlr3viz)
library(skimr)
library(DataExplorer)


###useful convert
# data = final1 %>% mutate_at(vars(one_of("host")), funs( as.factor)) 
# data = final1 %>% mutate_if(sapply(data_test, is.character), as.factor)
# dataset %>% mutate(across(where(~all(. %in% c(0, 1))), factor, labels = c(FALSE, TRUE))) %>% mutate_if(is.factor, as.logical) %>% mutate_at(vars("Res"), as.factor)

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

dt_select <- test_datasets %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(BMI = cut(BMI, breaks = c(-Inf, 28, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(HCTdALB = test_datasets$`@1d_hct`/test_datasets$`@1d_ALB`) %>%
   mutate(`@1d_ALB` = cut(`@1d_ALB`, breaks = c(-Inf, 25, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`24h血糖（mmol/L）` = cut(`24h血糖（mmol/L）`, breaks = c(-Inf, 14, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(Hypertension = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==1, "1", "0")) %>% 
   mutate(Diabetes = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==2, "1", "0")) %>% 
   mutate(Hypertension = case_when(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==3 ~ "1",TRUE ~ Hypertension)) %>% 
   mutate(Diabetes = case_when(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==3 ~ "1",TRUE ~Diabetes)) %>%
   mutate(Others_disease = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4其他)`==4, "1", "0")) %>% 
   select(-c("病案号","姓名","体重（kg","身高m","天数","中性粒细胞","吸入性损伤...46"))   %>% 
   rename_with(~str_remove(., '[@）]+')) %>% 
   mutate_if(is.character, as.factor) %>%  
   rename("耐受"="主要结局：耐受0，不耐受1") %>% rename("吸入性损伤"="吸入性损伤...12") 
   
colnames(dt_select) <- gsub('[(（].*','',colnames(dt_select))

dt_train = dt_select  %>% 
    select(c("耐受","年龄","性别","BMI","吸入性损伤","天数","24h血糖","TBSA","烧伤指数","III度","脓毒症","HCTdALB","1d_ALB","1dHB","1d_plt","1d_淋巴细胞","1d_TP","1d_TBIL","1d_DBIL","1d_CRE","1d_BUN","1d_PA","1d_LAC","Hypertension","Diabetes","Others_disease"))

colnames(dt_train)=c("Res","Age","Sex","Bmi","Inhalation_injury","Day","Glu_24h","TBSA","Burn_index","III_index","Sepsis","HCTdALB","ALB","HB","Plt","Lymphocyte","TP","TBIL","DBIL","CRE","BUN","PA","LAC","Hypertension","Diabetes","Others_disease")

validate_datasets <- read_excel("validate_datasets.xlsx")

dt_test = validate_datasets %>%
   mutate(年龄 = str_remove(年龄, "岁")) %>% mutate_at(vars("年龄"), as.numeric) %>%
   mutate(年龄 = cut(年龄, breaks = c(-Inf, 60, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(BMI = cut(BMI, breaks = c(-Inf, 28, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(HCTdALB = validate_datasets$`@1d_hct`/validate_datasets$`@1d_ALB`) %>% 
   mutate(`@1d_ALB` = cut(`@1d_ALB`, breaks = c(-Inf, 25, Inf), right = FALSE, labels = c("0", "1"))) %>%
   mutate(`24h血糖` = cut(`24h血糖`, breaks = c(-Inf, 14, Inf), right = FALSE, labels = c("0", "1"))) %>%
   # mutate(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`= ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==0, "0", "1")) %>%
   mutate(Hypertension = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==1, "1", "0")) %>% 
   mutate(Diabetes = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==2, "1", "0")) %>% 
   mutate(Hypertension = case_when(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==3 ~ "1",TRUE ~ Hypertension)) %>% 
   mutate(Diabetes = case_when(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==3 ~ "1",TRUE ~Diabetes)) %>%
   mutate(Others_disease = ifelse(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==5, "1", "0")) %>% 
   mutate(Others_disease = case_when(`基础疾病(0无1高血压2糖尿病3两者都有4肿瘤5其他)`==4 ~ "1",TRUE ~Others_disease))%>%
   mutate(天数 = cut(天数, breaks = c(-Inf, 4, Inf), right = FALSE, labels = c("0", "1"))) %>%
   select(c("主要结局：耐受0，不耐受1","年龄","性别","BMI","吸入性损伤","天数","24h血糖","TBSA","烧伤指数","III度","脓毒症（无0有1）","HCTdALB","@1d_ALB","@1dHB","@1d_plt","@1d_淋巴细胞","@1d_TP","@1d_TBIL","@1d_DBIL","@1d_CRE","@1d_BUN","@1d_PA","@1d_LAC","Hypertension","Diabetes","Others_disease")) %>% 
   mutate_if(is.character, as.factor) 
colnames(dt_test)=c("Res","Age","Sex","Bmi","Inhalation_injury","Day","Glu_24h","TBSA","Burn_index","III_index","Sepsis","HCTdALB","ALB","HB","Plt","Lymphocyte","TP","TBIL","DBIL","CRE","BUN","PA","LAC","Hypertension","Diabetes","Others_disease")

skimr::skim(dt_test)

dataset = rbind(dt_train,dt_test) # %>% select(c("Res", "Age","Inhalation_injury","Burn_index","HCTdALB","Plt","TBIL"))




###################data report separate


plot_str(dataset)

plot_missing(dataset)
# dataset <- drop_columns(dataset, "total missing")
plot_bar(dataset,nrow = 4L, ncol = 4L)
plot_histogram(dataset)
#select_if(negate(is.numeric))
qq_data <- dataset %>% select_if(is.numeric)
plot_qq(qq_data, sampled_rows = 1000L,nrow = 4L, ncol = 4L)

log_qq_data <- update_columns(qq_data, 1:ncol(qq_data), function(x) log(x + 1))

plot_qq(log_qq_data, sampled_rows = 1000L)

plot_qq(qq_data, by = "TBSA", sampled_rows = 1000L,nrow = 4L, ncol = 4L)

plot_correlation(na.omit(dataset), maxcat = 5L)

plot_boxplot(dataset, by = "Res",nrow = 4L, ncol = 4L)

###########data report 
create_report(dataset, y = "Res")


##############
raw_test = TaskClassif$new("raw_test", dataset, target = "Res")

factor_pipeline = po("encode", method = "treatment", affect_columns = selector_type("factor"))

lrn_xgb = lrn("classif.xgboost", nrounds = 100)

factor_pipeline$train(list(raw_test))[[1]]$data()







#########0,1 to logical
# dataset = dataset %>% mutate(across(where(~all(. %in% c(0, 1))), factor, labels = c(FALSE, TRUE))) %>% mutate_if(is.factor, as.logical) %>% mutate_at(vars("Res"), as.factor)

task = TaskClassif$new("shao_test", dataset, target = "Res")

factor_pipeline = po("encode", method = "treatment", affect_columns = selector_type("factor")) %>>% po("imputeoor") 
factor_pipeline$train(task)[[1]]$data()

task2 = factor_pipeline$train(task)[[1]]

###########
#A Resampling is instantiated for a task with a different number of observations
task2$set_col_roles(c("Age","Bmi"), "stratum")
table(task2$data(cols = c("Age","Bmi")))

rsmp_cv10 = rsmp("cv", folds = 10)
rsmp_cv10$instantiate(task2)

fold1 = prop.table(table(task2$data(rows = rsmp_cv10$test_set(1),cols = "Bmi")))
fold2 = prop.table(table(task2$data(rows = rsmp_cv10$test_set(2),cols = "Bmi")))
rbind("Fold 1" = fold1, "Fold 2" = fold2)


############
#custom = rsmp("custom")
#train_sets = list(1:nrow(dt_train), (nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)))
#test_sets = list((nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)), 1:nrow(dt_train))
#custom$instantiate(task2, train_sets, test_sets)

rsmp_tuner =  rsmp("cv", folds = 10)

custom1 = rsmp("custom")
train_sets = list(1:nrow(dt_train))
test_sets = list((nrow(dt_train)+1):(nrow(dt_train)+nrow(dt_test)))
custom1$instantiate(task2, train_sets, test_sets)

#task = TaskClassif$new("shao_test", dt_train, target = "Res")

at_featureless = auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.featureless", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_log_reg = auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.log_reg", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_cv_glmnet = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.cv_glmnet", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_lda= auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.lda", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_naive_bayes = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.naive_bayes", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_xgboost = auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.xgboost", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_ranger = auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.ranger", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_svm = auto_tuner(tuner=tnr("random_search"), learner =  lrn("classif.svm", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)
at_rpart = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.rpart", predict_type = "prob"),resampling = rsmp_tuner, measure = msr("classif.auc", average = "micro"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)

auto <- auto_fselector(
  fselector = fs("random_search"),
  learner = lrn("classif.log_reg", predict_type = "prob"),
  resampling = rsmp("loo"),
  measure = msr("classif.auc", average = "micro"),
  terminator = trm("evals", n_evals = 10)
)


learners <- c(auto,at_featureless,at_log_reg,at_cv_glmnet, at_lda,at_naive_bayes,at_xgboost,at_ranger,at_svm,at_rpart)
# learners = po("encode") %>>% learners
measures <- msrs(c("classif.auc", "classif.bacc", "classif.bbrier"))

#Benchmarking
set.seed(372)
design = benchmark_grid(tasks =task2, learners = learners, resamplings = custom1)
bmr = benchmark(design, store_models = TRUE)
bmr$aggregate(measures)

########
at_cv_glmnet$train(task2, train_sets[[1]])
prediction = at_cv_glmnet$predict(task2, test_sets[[1]])

prediction$confusion
prediction$score(msr("classif.auc"))

###########res unbalance
new_thresh = proportions(table(tsk_zoo$truth(splits$train)))
prediction$set_threshold(new_thresh)


#####weight
df = task2$data()
df$weights = ifelse(df$Age == 1, 2, 1)


# create new task and role
task2_weighted = as_task_classif(df, target = "Res")
task2_weighted$set_col_roles("weights", roles = "weight")

# compare weighted and unweighted predictions

prediction_weighted = at_cv_glmnet$train(task2_weighted, train_sets[[1]])
prediction_weighted = at_cv_glmnet$predict(task2, test_sets[[1]])
prediction_weighted$score(msr("classif.auc"))






results <- bmr$aggregate(measures, average = "micro")
print(results)
autoplot(bmr, measure = msr("classif.auc"))
autoplot(bmr, type = "roc")

########robust test
at_cv_glmnet_ro = auto_tuner(tuner=tnr("random_search"), learner = as_learner(ppl("robustify") %>>% lrn("classif.cv_glmnet", predict_type = "prob")),resampling = rsmp("cv", folds = 3), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)

design = benchmark_grid(tasks =task2, learners = c(at_cv_glmnet,at_cv_glmnet_ro), resamplings = rsmp("cv", folds = 5) )
bmr_new = benchmark(design, store_models = TRUE)
bmr$combine(bmr_new)

#########traning and test new data
at_cv_glmnet$train(task)

at_log_reg$predict_newdata(dataset)

####po stack run
library(mlr3verse)

task = TaskClassif$new("shao_test", dataset, target = "Res")

graph_stack = 
   po("scale") %>>%
   # po("filter", filter = flt("auc"), filter.frac = 0.5) %>>% # filter.frac应该做调参, 还有其它可选filter.nfeat, filter.cutoff, filter.permuted
   #gunion(list(
   #po("learner_cv", lrn("classif.cv_glmnet")),
   #po("learner_cv", lrn("classif.svm")),
   #po("nop"))) %>>%
   #po("featureunion") %>>%
   lrn("classif.ranger",predict_type = "prob")

learner = as_learner(graph_stack)
rr = resample(task, learner, rsmp("cv", folds = 5))

rr$score()
rr$aggregate(msrs("classif.auc"))


at_ranger = auto_tuner(tuner=tnr("random_search"), learner = lrn("classif.ranger", predict_type = "prob"),resampling = rsmp("cv", folds = 3), measure = msr("classif.auc"),term_evals = 20,store_tuning_instance = TRUE,store_models = TRUE)


design = benchmark_grid(task, list(at_ranger, as_learner(graph_stack)),
  rsmp("cv", folds = 5))
bmr = benchmark(design, store_models = TRUE)
bmr$aggregate()[, .(learner_id, classif.ce)]
bmr$aggregate(msr("classif.auc"))[, .(learner_id, classif.auc)]

#####################

# create mlr h2o model
library(mlr)

dat=rbind(dt_train,dt_test)

train=1:nrow(dt_train)

datTrain <- dat[train, ]
datTest <- dat[-train, ]

task2 <- makeClassifTask(data =  datTrain , target = "Res")

learner <- makeLearner("classif.h2o.deeplearning", predict.type = "prob", 
                       par.vals = list(reproducible = TRUE,
                                       seed = 1))
Mod <- train(learner, task2)

# Test predictions
pred <- predict(Mod, newdata = datTest)
# Evaluate performance accuracy & area under curve 
performance(pred, measures = list(acc, auc)) 

set.seed(1234)
# Tune epoch parameter
param_set <- makeParamSet(
  makeNumericParam("epochs", lower = 1, upper = 10))
rdesc <- makeResampleDesc("CV", iters = 3L, predict = "both") 

ctrl <- makeTuneControlRandom(maxit = 3)

res <- tuneParams(
  learner = learner, task = task2, resampling = rdesc, measures = list(auc, acc),
  par.set = param_set, control = ctrl
)

resample(learner, task2, cv3, list(auc, acc))

set.seed(1234)
# plugging the tuned value into model and checking performance again:
learner <- makeLearner("classif.h2o.deeplearning", predict.type = "prob", 
                       par.vals = list(epochs = 4.54,
                                       reproducible = TRUE,
                                       seed = 1))
Mod <- train(learner, task2)

# Test predictions
pred1 <- predict(Mod, newdata = datTest)
# Evaluate performance accuracy & area under curve 
performance(pred1, measures = list(acc, auc))

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