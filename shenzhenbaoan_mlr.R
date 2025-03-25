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
library(easystats)

library(readr)
library(dplyr)
library(glmnet)
library(mlr3verse)
library(iml)

all_datasets <- read_csv("26_genes.csv")

dataset <- all_datasets %>% filter(Group!="LTBI") %>% select(-c("THRSP")) %>%
            rename(Res="Group")  %>% mutate_at(vars("Res"), as.factor)

# create_report(dataset, y = "Res")
tsk = TaskClassif$new("baoan", dataset, target = "Res")

library(mlr3filters)
flt("information_gain")$calculate(tsk) %>% as.data.table() %>% write_csv("importance.csv")
flt("correlation", method = "kendall")$calculate(tsk) %>% as.data.table() %>% write_csv("correlation.csv")
# flts(c("mrmr", "jmim"))$calculate(tsk) %>% as.data.table() %>% write_csv("mrmr+jmim.csv")
flt("mrmr")$calculate(tsk) %>% as.data.table() %>% write_csv("mrmr.csv")
##########


split = partition(tsk)
learner.lasso.cv <- lrn("classif.cv_glmnet")

# lasso_gr = po('encode') %>>% po('scale') %>>% po(learner.lasso.cv)
lasso_gr = po("imputeoor") %>>% po(learner.lasso.cv)
lasso_glrn = GraphLearner$new(lasso_gr)


lasso_glrn$train(tsk, row_ids = split$train)

credit_x = tsk$data(rows = split$test,
  cols = tsk$feature_names)
# target in test data
credit_y = tsk$data(rows = split$test,
  cols = tsk$target_names)

predictor = Predictor$new(lasso_glrn, data = credit_x, y = credit_y)


importance = FeatureImp$new(predictor, loss = "ce", n.repetitions = 10)
importance$plot()




###############

lrn = lrn("classif.rpart", predict_type = "prob")
fit = train(lrn, tsk)
imp = generateFeatureImportanceData(tsk, "permutation.importance",
  lrn, "Petal.Width", nmc = 10L, local = TRUE)









###################data report separate

instance = fselect(
  # fselector = fs("sequential"),
  fselector = fs("random_search", batch_size = 5),
  task = TaskClassif$new("baoan", dataset, target = "Res"),
  learner = lrn("classif.rpart"),
  # learner =  po("imputeoor") %>>% lrn("classif.cv_glmnet"),
  resampling = rsmp("cv", folds = 10),
  # resampling = rsmp("holdout"),
  measures = msr("classif.ce")
  # terminator = trm("evals", n_evals = 20)
)

instance$result

fselector$optimize(instance)

instance$result_feature_set


task = tsk("pima")

# load learner
learner = lrn("classif.rpart")

# feature selection on the pima indians diabetes data set
instance = fselect(
  fselector = fs("sequential"),
  task = task,
  learner = learner,
  resampling = rsmp("holdout"),
  measure = msr("classif.ce")
)

# best performing feature subset
instance$result