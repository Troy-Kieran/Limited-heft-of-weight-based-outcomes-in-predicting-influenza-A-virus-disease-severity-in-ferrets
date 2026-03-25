########################################################################################
########################################################################################
### Morbidity2 Analysis - Machine Learning & Other Models
### 10 December 2025 - 15 January 2026
### Troy J. Kieran 
########################################################################################
########################################################################################
### Load packages

library(tidyverse)
library(tidylog)
library(caret)
library(gbm)
library(funModeling)

########################################################################################
### Import Data
fullMorbidity <- ## data generated from 1_morbidity2_data_cleaning.R

########################################################################################

## standing function for machine learning
fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## function for Matthew's Correlation Coeffecient
mcc_func <- function(TP, FP, FN, TN){
  ((TP*TN) - (FP*FN))/
    sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))}

########################################################################################

## subset to relevant variables to test
vars <- fullMorbidity %>% 
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f, 
                  temp_5, wt_loss_slope_max, wt_rmssd_rate, wt_acf1, 
                  temp_rmssd_rate, temp_acf1)) %>% 
  drop_na()

########################################################################################

## categorize new weight loss metrics into high vs not

## check initial bins
describe(equal_freq(vars$wt_loss_slope_max, n_bins = 3))
describe(equal_freq(vars$wt_rmssd_rate, n_bins = 3))
describe(equal_freq(vars$wt_acf1, n_bins = 3))

vars2 <- vars %>%
  mutate(wt_loss_slope_max_cat = case_when(wt_loss_slope_max >= -1.05 ~ 'low',
                                           wt_loss_slope_max >= -2.08 ~ 'med',
                                           wt_loss_slope_max > -6.4 ~ 'high'),
         wt_rmssd_rate_cat = case_when(wt_rmssd_rate < 1.709 ~ 'low',
                                       wt_rmssd_rate <= 2.31 ~ 'med',
                                       wt_rmssd_rate <= 6.61 ~ 'high'),
         wt_acf1_cat = case_when(wt_acf1 <= 0.673 ~ 'low',
                                 wt_acf1 <= 0.882 ~ 'med',
                                 wt_acf1 <= 1 ~ 'high')) %>%
  mutate(wt_loss_high_slope = ifelse(wt_loss_slope_max_cat == 'high', 'yes', 'no'),
         wt_loss_high_rmssd = ifelse(wt_rmssd_rate_cat == 'high', 'yes', 'no'),
         wt_loss_high_acf1 = ifelse(wt_acf1_cat == 'high', 'yes', 'no')) %>%
  select(-contains('_cat'))


########################################################################################

## Elastic Net model
OUT_COL <- 'wt_loss_high'

X1 <- subset(vars, select = setdiff(names(vars), OUT_COL))
## convert all character columns to factors 
X1[] <- lapply(X1, function(x) if(is.character(x)) as.factor(x) else x)
## identify factor variables with at least 2 levels
factor_vars <- sapply(X1, function(x) is.factor(x) && nlevels(x) >= 2)

XX <- model.matrix(~ ., data = X1)
yy <- vars[[OUT_COL]]

set.seed(9595)
elastic_model <- glmnet::cv.glmnet(x = XX, y = yy, alpha = 0.35, nfolds = 10, 
                                   family = "binomial")
best_lambda <- elastic_model$lambda.min

best_model <- glmnet::glmnet(XX, y = yy, alpha = 0.35, lambda = best_lambda, 
                             family = "binomial")

all_coefs <- capture.output(coef(best_model))
data.frame(Coef = all_coefs)

#####################################

## rotate through outcome
OUT_COL <- 'wt_loss_high_slope'
OUT_COL <- 'wt_loss_high_rmssd'
OUT_COL <- 'wt_loss_high_acf1'

X1 <- subset(vars2, select = setdiff(names(vars2), OUT_COL))
## convert all character columns to factors 
X1[] <- lapply(X1, function(x) if(is.character(x)) as.factor(x) else x)
## identify factor variables with at least 2 levels
factor_vars <- sapply(X1, function(x) is.factor(x) && nlevels(x) >= 2)

XX <- model.matrix(~ ., data = X1)
yy <- vars2[[OUT_COL]]

set.seed(9595)
elastic_model <- glmnet::cv.glmnet(x = XX, y = yy, alpha = 0.35, nfolds = 10, 
                                   family = "binomial")
best_lambda <- elastic_model$lambda.min

best_model <- glmnet::glmnet(XX, y = yy, alpha = 0.35, lambda = best_lambda, 
                             family = "binomial")

all_coefs <- capture.output(coef(best_model))
data.frame(Coef = all_coefs)

########################################################################################

## select variables for each test

base <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                          temp_5))

## base model + new weight metrics
test1 <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                          temp_5, wt_loss_slope_max))

test2 <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_5, wt_rmssd_rate))

test3 <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_5, wt_acf1))

## base but swap out temp_5 for other temp metrics
test4 <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_rmssd_rate))

test5 <- vars %>% select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_acf1))


## swap out weight related outcome variable
base1 <- vars2 %>% select(c(wt_loss_high_slope, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                          temp_5))

base2 <- vars2 %>% select(c(wt_loss_high_rmssd, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_5))

base3 <- vars2 %>% select(c(wt_loss_high_acf1, HPAI_MBAA, RBS, PBS, HA, AUC_6_f,
                           temp_5))

########################################################################################

### prep data for machine learning

## dummy variables
base_dummy <- fastDummies::dummy_cols(
  base, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

test1_dummy <- fastDummies::dummy_cols(
  test1, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

test2_dummy <- fastDummies::dummy_cols(
  test2, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

test3_dummy <- fastDummies::dummy_cols(
  test3, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

test4_dummy <- fastDummies::dummy_cols(
  test4, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

test5_dummy <- fastDummies::dummy_cols(
  test5, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

base1_dummy <- fastDummies::dummy_cols(
  base1, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

base2_dummy <- fastDummies::dummy_cols(
  base2, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

base3_dummy <- fastDummies::dummy_cols(
  base3, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

## outcome a factor
base_dummy$wt_loss_high <- as.factor(base_dummy$wt_loss_high)
test1_dummy$wt_loss_high <- as.factor(test1_dummy$wt_loss_high)
test2_dummy$wt_loss_high <- as.factor(test2_dummy$wt_loss_high)
test3_dummy$wt_loss_high <- as.factor(test3_dummy$wt_loss_high)
test4_dummy$wt_loss_high <- as.factor(test4_dummy$wt_loss_high)
test5_dummy$wt_loss_high <- as.factor(test5_dummy$wt_loss_high)
base1_dummy$wt_loss_high_slope <- as.factor(base1_dummy$wt_loss_high_slope)
base2_dummy$wt_loss_high_rmssd <- as.factor(base2_dummy$wt_loss_high_rmssd)
base3_dummy$wt_loss_high_acf1 <- as.factor(base3_dummy$wt_loss_high_acf1)

## train/test data split
set.seed(9595)
base_dumSplit <- rsample::initial_split(base_dummy, prop = 0.70)
test1_dumSplit <- rsample::initial_split(test1_dummy, prop = 0.70)
test2_dumSplit <- rsample::initial_split(test2_dummy, prop = 0.70)
test3_dumSplit <- rsample::initial_split(test3_dummy, prop = 0.70)
test4_dumSplit <- rsample::initial_split(test4_dummy, prop = 0.70)
test5_dumSplit <- rsample::initial_split(test5_dummy, prop = 0.70)
base1_dumSplit <- rsample::initial_split(base1_dummy, prop = 0.70)
base2_dumSplit <- rsample::initial_split(base2_dummy, prop = 0.70)
base3_dumSplit <- rsample::initial_split(base3_dummy, prop = 0.70)

trainData_base <- rsample::training(base_dumSplit)
testData_base <- rsample::testing(base_dumSplit)
trainData_test1 <- rsample::training(test1_dumSplit)
testData_test1 <- rsample::testing(test1_dumSplit)
trainData_test2 <- rsample::training(test2_dumSplit)
testData_test2 <- rsample::testing(test2_dumSplit)
trainData_test3 <- rsample::training(test3_dumSplit)
testData_test3 <- rsample::testing(test3_dumSplit)
trainData_test4 <- rsample::training(test4_dumSplit)
testData_test4 <- rsample::testing(test4_dumSplit)
trainData_test5 <- rsample::training(test5_dumSplit)
testData_test5 <- rsample::testing(test5_dumSplit)
trainData_base1 <- rsample::training(base1_dumSplit)
testData_base1 <- rsample::testing(base1_dumSplit)
trainData_base2 <- rsample::training(base2_dumSplit)
testData_base2 <- rsample::testing(base2_dumSplit)
trainData_base3 <- rsample::training(base3_dumSplit)
testData_base3 <- rsample::testing(base3_dumSplit)

########################################################################################

### train machine learning

set.seed(2626)
mbase <- train(wt_loss_high ~ ., data = trainData_base,
               method = "gbm",
               metric = "Balanced_Accuracy",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
mtest1 <- train(wt_loss_high ~ ., data = trainData_test1,
               method = "gbm",
               metric = "Balanced_Accuracy",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
mtest2 <- train(wt_loss_high ~ ., data = trainData_test2,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

set.seed(2626)
mtest3 <- train(wt_loss_high ~ ., data = trainData_test3,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

set.seed(2626)
mtest4 <- train(wt_loss_high ~ ., data = trainData_test4,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

set.seed(2626)
mtest5 <- train(wt_loss_high ~ ., data = trainData_test5,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

set.seed(2626)
mbase1 <- train(wt_loss_high_slope ~ ., data = trainData_base1,
               method = "gbm",
               metric = "Balanced_Accuracy",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
mbase2 <- train(wt_loss_high_rmssd ~ ., data = trainData_base2,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

set.seed(2626)
mbase3 <- train(wt_loss_high_acf1 ~ ., data = trainData_base3,
                method = "gbm",
                metric = "Balanced_Accuracy",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

###

## resample metrics
resamps <- resamples(list(Base = mbase,
                          wt_loss_slope_max = mtest1,
                          wt_rmssd_rate = mtest2,
                          wt_acf1 = mtest3,
                          temp_rmssd_rate = mtest4,
                          temp_acf1 = mtest5))

resamps <- resamples(list(Base = mbase,
                          wt_loss_slope_max = mbase1,
                          wt_rmssd_rate = mbase2,
                          wt_acf1 = mbase3))

summary(resamps)
bwplot(resamps)
dotplot(resamps)
modelCor(resamps)
splom(resamps)

## check for variable importance
#ggplot(varImp(mbase))
varImp(mbase)
varImp(mtest1)
varImp(mtest2)
varImp(mtest3)
varImp(mtest4)
varImp(mtest5)

varImp(mbase1)
varImp(mbase2)
varImp(mbase3)

########################################################################################

### test machine learning

## rotate through
testData <- testData_base
model <- mbase
#testData <- testData_test1
#model <- mtest1
#testData <- testData_test2
#model <- mtest2
#testData <- testData_test3
#model <- mtest3
#testData <- testData_test4
#model <- mtest4
#testData <- testData_test5
#model <- mtest5

###

## check a range of probability values
thresholds <- seq(0.5, 0.95, by = 0.05)

## empty list to populate
prob_results <- list()
## function
for (threshold in thresholds) {
  probTest <- predict(model, testData, type = 'prob')
  probTest <- factor(ifelse(probTest$no >= threshold, 'no', 'yes')) %>%
    as.data.frame()
  
  probTruth <- testData %>%
    na.omit() %>%
    dplyr::select(wt_loss_high) %>%
    cbind(., probTest)
  
  confusion_matrix <- confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")
  
  prob_results[[as.character(threshold)]] <- list(
    threshold = threshold,
    confusion_matrix = confusion_matrix)}

## access the results for each threshold, e.g., results[['0.5']]
prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]


## Matthew's Correlation Coefficient
## base - 0.55 prob
mcc_func(157, 25, 21, 29) ## 0.4306871
## test1 - 0.65 prob - wt_loss_slope_max
mcc_func(144, 10, 41, 37) ## 0.4812227
## test2 - 0.85 prob - wt_rmssd_rate
mcc_func(120, 13, 55, 44) ## 0.3983281
## test3 - 0.70 prob - wt_acf1
mcc_func(149, 11, 24, 48) ## 0.6352108
## test4 - 0.75 prob - temp_rmssd_rate
mcc_func(128, 21, 46, 37) ## 0.3374586
## test5 - 0.75 prob - temp_acf1
mcc_func(144, 31, 32, 25) ## 0.2630273

########################################################################################

## rotate through
#testData <- testData_base1
#model <- mbase1
testData <- testData_base2
model <- mbase2
#testData <- testData_base3
#model <- mbase3

###

## check a range of probability values
thresholds <- seq(0.5, 0.95, by = 0.05)

## empty list to populate
prob_results <- list()
## function - swap out outcome variable for relevant one
for (threshold in thresholds) {
  probTest <- predict(model, testData, type = 'prob')
  probTest <- factor(ifelse(probTest$no >= threshold, 'no', 'yes')) %>%
    as.data.frame()
  
  probTruth <- testData %>%
    na.omit() %>%
    dplyr::select(wt_loss_high_rmssd) %>%
    cbind(., probTest)
  
  confusion_matrix <- confusionMatrix(probTruth$., probTruth$wt_loss_high_rmssd, mode = "everything")
  
  prob_results[[as.character(threshold)]] <- list(
    threshold = threshold,
    confusion_matrix = confusion_matrix)}

## access the results for each threshold, e.g., results[['0.5']]
prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]


## Matthew's Correlation Coefficient
## base - 0.55 prob - wt_loss_high
mcc_func(157, 25, 21, 29) ## 0.4306871
## base1 - 0.60 prob - wt_loss_high_slope
mcc_func(122, 37, 35, 38) ## 0.2857828
## base2 - 0.60 prob - wt_loss_high_rmssd
mcc_func(119, 41, 40, 32) ## 0.1874878
## base3 - 0.60 prob - wt_loss_high_acf1
mcc_func(118, 51, 26, 37) ## 0.2617203

########################################################################################
########################################################################################

### Explore weight/temp/titer associations

## Fill in the d2,4,6,8 inoc with the mean of day on either side
## Convert NAs to median by units+Virus for dx_inoc
## Convert dx_inoc values back to NA for euth individuals
fullMorbidity_nw_imputed <- 
  fullMorbidity %>%
  dplyr::filter(samp_type == "NW") %>%
  drop_na(wt_0) %>%
  mutate(d2_inoc = ifelse(is.na(d2_inoc), 
                          rowMeans(.[, c("d1_inoc", "d3_inoc")]), 
                          d2_inoc),
         d4_inoc = ifelse(is.na(d4_inoc), 
                          rowMeans(.[, c("d3_inoc", "d5_inoc")]), 
                          d4_inoc),
         d6_inoc = ifelse(is.na(d6_inoc), 
                          rowMeans(.[, c("d5_inoc", "d7_inoc")]), 
                          d6_inoc),
         d8_inoc = ifelse(is.na(d8_inoc), 
                          rowMeans(.[, c("d7_inoc", "d9_inoc")]), 
                          d8_inoc),
         d3_inoc = ifelse(is.na(d3_inoc), 
                          rowMeans(.[, c("d2_inoc", "d4_inoc")]), 
                          d3_inoc),
         d5_inoc = ifelse(is.na(d5_inoc), 
                          rowMeans(.[, c("d4_inoc", "d6_inoc")]), 
                          d5_inoc),
         d7_inoc = ifelse(is.na(d7_inoc), 
                          rowMeans(.[, c("d6_inoc", "d8_inoc")]), 
                          d7_inoc)) %>%
  group_by(units, Virus) %>% 
  mutate(d1_inoc = ifelse(is.na(d1_inoc),
                          median(d1_inoc, na.rm = TRUE),
                          d1_inoc),
         d2_inoc = ifelse(is.na(d2_inoc),
                          median(d2_inoc, na.rm = TRUE),
                          d2_inoc),
         d3_inoc = ifelse(is.na(d3_inoc),
                          median(d3_inoc, na.rm = TRUE),
                          d3_inoc),
         d4_inoc = ifelse(is.na(d4_inoc),
                          median(d4_inoc, na.rm = TRUE),
                          d4_inoc),
         d5_inoc = ifelse(is.na(d5_inoc),
                          median(d5_inoc, na.rm = TRUE),
                          d5_inoc),
         d6_inoc = ifelse(is.na(d6_inoc),
                          median(d6_inoc, na.rm = TRUE),
                          d6_inoc),
         d7_inoc = ifelse(is.na(d7_inoc),
                          median(d7_inoc, na.rm = TRUE),
                          d7_inoc),
         d8_inoc = ifelse(is.na(d8_inoc),
                          median(d8_inoc, na.rm = TRUE),
                          d8_inoc),
         d9_inoc = ifelse(is.na(d9_inoc),
                          median(d9_inoc, na.rm = TRUE),
                          d9_inoc)) %>% 
  ungroup() %>% 
  mutate(d5_inoc = ifelse(d5_inoc_euth == "yes", NA, d5_inoc),
         d6_inoc = ifelse(d5_inoc_euth == "yes", NA, d6_inoc),
         d7_inoc = ifelse(d7_inoc_euth == "yes", NA, d7_inoc),
         d8_inoc = ifelse(d7_inoc_euth == "yes", NA, d8_inoc),
         d9_inoc = ifelse(d9_inoc_euth == "yes", NA, d9_inoc))

###

## trim some unneeded columns
fullMorbidity_nw_imputed <- fullMorbidity_nw_imputed %>%
  select(-contains(c('inoc_euth', 'DC_euth', 'RD_euth', '_rep', 
                     '_avg', 'AUC', '_titer', '_sero', '_LOD', 
                     'peak', 'slope', '_max', 'NW_typical', 'samp_type'))) %>% 
  select(where(~ !all(is.na(.))))

fullMorbidity_nw_imputed2 <- fullMorbidity_nw_imputed %>%
  select(-contains(c('inoc_euth', 'DC_euth', 'RD_euth', '_rep', 
                     '_avg', 'AUC', '_titer', '_sero', '_LOD', 
                     'peak', 'NW_typical', 'samp_type'))) %>% 
  select(where(~ !all(is.na(.))))

###

df <- fullMorbidity_nw_imputed
id_col <- "Ferret"

## weight long
wt_long <- df %>%
  pivot_longer(
    cols = matches("^wt_\\d+(_n)?$"),
    names_to = "name",
    values_to = "wt") %>%
  mutate(
    day = as.integer(str_extract(name, "\\d+")),
    wt_norm = str_detect(name, "_n")) %>%
  select(!!id_col, Virus, day, wt, wt_norm)

## temperature long
tp_long <- df %>%
  pivot_longer(
    cols = matches("^tp_\\d+(_n)?$"),
    names_to = "name",
    values_to = "tp") %>%
  mutate(
    day = as.integer(str_extract(name, "\\d+")),
    tp_norm = str_detect(name, "_n")) %>%
  select(!!id_col, Virus, day, tp, tp_norm)

## titers long (days embedded as d#_inoc)
titer_long <- df %>%
  pivot_longer(
    cols = matches("^d\\d+_inoc$"),
    names_to = "name",
    values_to = "titer") %>%
  mutate(
    day = as.integer(str_extract(name, "\\d+"))) %>%
  select(!!id_col, Virus, day, titer)

## outcome and other subject-level vars
subject_vars <- df %>%
  select(all_of(c(id_col, "Virus", "lethal"))) %>%
  distinct()

long <- wt_long %>%
  full_join(tp_long, by = c(id_col, "Virus", "day")) %>%
  full_join(titer_long, by = c(id_col, "Virus", "day")) %>%
  left_join(subject_vars, by = c(id_col, "Virus"))

long <- long %>%
  ## keep raw wt/tp
  #filter(!wt_norm %in% TRUE | is.na(wt_norm)) %>%
  #filter(!tp_norm %in% TRUE | is.na(tp_norm)) %>%
  ## keep normalized wt/tp
  filter(wt_norm %in% TRUE | is.na(wt_norm)) %>%
  filter(tp_norm %in% TRUE | is.na(tp_norm)) %>%
  select(-wt_norm, -tp_norm) %>%
  arrange(.data[[id_col]], day) %>% 
  rename(weight = wt, temp = tp)

###

## distribution summaries by Virus and day
long %>%
  group_by(Virus, day) %>%
  dplyr::summarize(
    n = n(),
    mean_weight = mean(weight, na.rm = TRUE),
    sd_weight = sd(weight, na.rm = TRUE),
    mean_temp = mean(temp, na.rm = TRUE),
    sd_temp = sd(temp, na.rm = TRUE),
    mean_titer = mean(titer, na.rm = TRUE),
    sd_titer = sd(titer, na.rm = TRUE)) %>%
  arrange(Virus, day)

## repeated-measures correlation
rmcorr::rmcorr(Ferret, weight, titer, long) 
## raw: r = 0.3824007,  p = 7.493986e-211; norm: r = 0.3774142,  p = 5.266994e-205
rmcorr::rmcorr(Ferret, weight, temp, long)  
## raw: r = 0.01957248, p = 0.05647915;    norm: r = 0.01266895, p = 0.2170133
rmcorr::rmcorr(Ferret, titer, temp, long)   
## raw: r = 0.2310591,  p = 1.866462e-73;  norm: r = 0.2310591,  p = 1.866462e-73

########################################################################################

### longitudinal modeling: linear mixed-effects models
library(lme4)
library(lmerTest) 

## weight model: day effects, Virus differences & ferret within virus diff, and association with temp/titer
## using normalize weight/temp
m_weight <- lmer(
  weight ~ day + temp + titer + (day | Virus/Ferret),
  data = long, control = lmerControl(optimizer = "bobyqa"))
summary(m_weight)
report::report(m_weight)
# day = –0.97 (p < 2e‑16)
#   weight decreases by ~1 unit per standardized day.
#   very strong, highly significant decline over time.
# temp = +0.109 (p = 0.0054)
#   higher temperature is associated with slightly higher weight.
#   effect is small but statistically significant.
# titer = –0.204 (p = 6.5e‑14)
#   higher viral titer predicts lower weight.
#   strong negative association, highly significant.
# Virus-level variation:
#   Intercept SD = 1.36, Slope SD = 0.88
#   viruses differ modestly in baseline weight and in rate of weight change.
# Ferret-within-virus variation:
#   Intercept SD = 1.94, Slope SD = 0.67
#   individual ferrets differ more than viruses do—especially in baseline weight.
# Overall:
#   weight declines over time, declines with increasing viral load, and increases modestly with temperature.
#   individual variation bigger than viral differences.

## temp model
m_temp <- lmer(
  temp ~ day + weight + titer + (day | Virus/Ferret),
  data = long, control = lmerControl(optimizer = "bobyqa"))
summary(m_temp)
report::report(m_temp)
# day = –0.089 (p = 1.6e‑12)
#   temperature decreases over time.
#   effect is small but highly significant.
#   suggests a consistent downward trend in temperature as infection progresses.
# weight = +0.010 (p = 0.00059)
#   higher weight is associated with slightly higher temperature.
#   effect is modest but statistically significant.
# titer = –0.0025 (p = 0.77)
#   titer is not significantly associated with temperature after accounting for day and weight.
#   suggests temperature is not strongly driven by viral load.
# Virus-level variation:
#   Intercept SD = 0.56, Slope SD = 0.10
#   viruses differ moderately in baseline temperature and slightly in how temperature changes over time.
# Ferret-within-virus variation:
#   Intercept SD = 0.45, Slope SD = 0.035
#   individual ferrets differ in baseline temperature, but their day‑to‑day temperature slopes vary only minimally.
# Overall:  
#   temperature declines over time, is modestly higher in heavier animals, and shows no meaningful association with viral titer.
#   virus variation is moderate, ferret variation minimal.

## titer model (only days 1–9)
m_titer <- lmer(
  titer ~ day + weight + temp + (day | Virus/Ferret),
  control = lmerControl(optimizer = "bobyqa"),
  data = long %>% filter(day <= 9))
summary(m_titer)
report::report(m_titer)
# day = –0.592 (p < 2e‑16)
#   titer declines sharply over time.
#   very strong effect, consistent with viral clearance after peak infection.
#   magnitude is much larger than in the temperature model, indicating titer is highly time‑dependent.
# weight = –0.023 (p = 8.8e‑15)
#   higher weight is associated with lower titer.
#   robust negative relationship, consistent with the idea that animals losing more weight have higher viral burden.
# temp = +0.010 (p = 0.50)
#     temperature is not significantly associated with titer after accounting for day and weight.
#     mirrors the temperature model, where titer also failed to predict temperature.
# Virus-level variation:
#   Intercept SD = 1.28, Slope SD = 0.19
#   viruses differ substantially in baseline titer and moderately in how titer changes over time.
#   Correlation = –0.85
#   strong negative correlation suggests viruses with higher baseline titers tend to show faster declines.
# Ferret-within-virus variation:
#   Intercept SD = 0.41, Slope SD = 0.034
#   individuals within a virus tend to have similar trajectories. 
# Overall:  
#   titer drops strongly over time, is higher in animals with lower weight, and is not meaningfully related to temperature.
#   substantial virus variation, essentially no individual ferret variation.

## diagnostics for mixed models
simulationOutput <- DHARMa::simulateResiduals(m_weight)
simulationOutput <- DHARMa::simulateResiduals(m_temp)
simulationOutput <- DHARMa::simulateResiduals(m_titer)
plot(simulationOutput)

###

df2 <- fullMorbidity %>% drop_na(NW_typical) %>% 
  select(Ferret, Virus, AUC_6_f, wt_loss_slope_max, 
         wt_rmssd_rate, wt_acf1, temp_rmssd_rate, temp_acf1) %>% 
  filter(!if_all(wt_loss_slope_max:wt_acf1, is.na))


m_weight2 <- lmer(
  wt_rmssd_rate ~ temp_rmssd_rate + AUC_6_f + (1 | Virus),
  data = df2, control = lmerControl(optimizer = "bobyqa"))
summary(m_weight2)
report::report(m_weight2)
# temp_rmssd_rate = +0.6213 (p = 5.54e-12)
#   higher temp variability is associated with modestly higher weight variability.
#   statistically significant.
# AUC_6_f = 5.806e-03 (p = 0.487)
#   titer no effect, not significant
# Virus-level variation:
#   Intercept SD = 0.4534
#   viruses differ somewhat in weight variability.
# Overall:
#   weight and temp variability associated with no titer effect, modest virus variability

m_temp2 <- lmer(
  temp_rmssd_rate ~ wt_rmssd_rate + AUC_6_f + (1 | Virus),
  data = df2, control = lmerControl(optimizer = "bobyqa"))
summary(m_temp2)
report::report(m_temp2)
# wt_rmssd_rate = +0.091542  (p = 1.73e-11)
#   higher weight variability is associated with minor higher temp variability.
#   statistically significant.
# AUC_6_f = -0.002593  (p = 0.428)
#   titer no effect, not significant
# Virus-level variation:
#   Intercept SD = 0.1869 
#   viruses differ slightly in temp variability.
# Overall:
#   weight and temp variability associated with no titer effect, minor virus variability

m_titer2 <- lmer(
  AUC_6_f ~ wt_rmssd_rate + temp_rmssd_rate + (1 | Virus),
  control = lmerControl(optimizer = "bobyqa"),
  data = df2)
summary(m_titer2)
report::report(m_titer2)
# temp_rmssd_rate = -0.237491 (p = 0.506 )
#   higher temp variability is slightly negatively associated with higher titer.
#   not significant.
# wt_rmssd_rate = -0.002958 (p = 0.983 )
#   titer no effect, not significant
# Virus-level variation:
#   Intercept SD = 3.460
#   viruses differ greatly in titer AUC.
# Overall:
#   weight and temp variability not associated with titer AUC, large virus variability

########################################################################################

## summarize models for output
modelsummary::modelsummary(
  list("weight" = m_weight,
       "temp" = m_temp,
       "titer" = m_titer,
       "weight rmssd" = m_weight2,
       "temp rmssd" = m_temp2,
       "titer AUC 6" = m_titer2), 
  output = "markdown")

########################################################################################

## plot trends

long %>%
  filter(weight < 20) %>%
  ggplot(aes(x = weight, y = titer, color = day)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_viridis_c(option = "plasma", direction = -1, end = 0.9) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = "Relationship Between Weight and Viral Titer",
       x = "Weight", y = "Viral Titer", color = "Day")

long %>%
  filter(temp > -5) %>%
  ggplot(aes(x = weight, y = temp, color = day)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_viridis_c(option = "plasma", direction = -1, end = 0.9) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = "Relationship Between Weight and Temperature",
       x = "Weight", y = "Temperature", color = "Day")

long %>%
  filter(weight > -5) %>%
  ggplot(aes(x = temp, y = titer, color = day)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_viridis_c(option = "plasma", direction = -1, end = 0.9) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = "Relationship Between Temperature and Viral Titer",
       x = "Temperature", y = "Viral Titer", color = "Day")

########################################################################################
########################################################################################

### Explore Coefficient of Variation of Variables

## function to summarize variables for variation
summarize_variability <- function(data, var) {
  data %>%
    filter(NW_typical == "yes") %>%
    group_by(wt_loss_high, units) %>%
    reframe(
      med = median({{var}}, na.rm = TRUE),
      mean = mean({{var}}, na.rm = TRUE),
      sd = sd({{var}}, na.rm = TRUE),
      Q1 = quantile({{var}}, 0.25, na.rm = TRUE),
      Q3 = quantile({{var}}, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      mad_val = mad({{var}}, constant = 1, na.rm = TRUE),
      ### measures of data variability:
      ## Gini coefficient - measure of variability, 0-1
      gini = DescTools::Gini({{var}}, na.rm = TRUE),
      ## Coefficient of Variation (CV)
      cv = (sd / mean) * 100,
      ## Robust CV based on IQR (RCV𝑄)
      rcvq = (0.75 * (IQR/med)) * 100,
      ## Robust CV based on MAD (RCV𝑀)
      rcvm = (1.4826 * (mad_val/med)) * 100) %>%
    pivot_longer(cols = c(cv, rcvq, rcvm),
                 names_to = "metric", values_to = "value") %>% 
    drop_na(value) #%>%
    #filter(inoc_dose %in% c(5, 6, 7)) 
    }

## use function to plot, swapping out var
summarize_variability(data = fullMorbidity, 
                      var = wt_loss) %>%
  ggplot(aes(x = wt_loss_high, y = value, fill = metric)) +
  ## Amacero colors from the poisonfrogs R package - Peru
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
           fill = "#C6E74B", alpha = 0.3) +
  geom_hline(yintercept = 25, linetype = 'dashed',linewidth = 1.2, color = '#C6E74B') +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(labels = c(cv = "Coefficient of Variation (CV)", 
                               rcvq = "CV via IQR", rcvm = "CV via MAD"),
                    values = c("cv" = "#3981BE", "rcvq" = "#64B936", "rcvm" = "#7D2A1D")) +
  facet_grid(~ units, scales = "free_x", space = "free_x",
             labeller = labeller(units = c(EID = "Eggs", London = "Cells"),
                                 inoc_dose = c('5' = '5', '6' = '6', '7' = '7'))) +
  labs(x = "Subtype", y = "Coefficient of Variation", fill = element_blank()) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom',
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = alpha("#C6E74B", 0.3), color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))

###

## summarize one CV metric for multiple variables
summarize_rcvm <- function(data, vars) {
  data %>%
    filter(NW_typical == "yes") %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "val") %>%
    group_by(wt_loss_high, variable) %>%
    reframe(
      med = median(val, na.rm = TRUE),
      mad_val = mad(val, constant = 1, na.rm = TRUE),
      rcvm = ifelse(med == 0 | is.na(med), NA_real_, abs(1.4826 * (mad_val / med) * 100))) %>% drop_na()}

#summarize_rcvm(fullMorbidity, vars = c("wt_loss_slope_max", "wt_rmssd_rate", "wt_acf1"))

## use function to plot, swapping out multiple vars
summarize_rcvm(data = fullMorbidity, 
               vars = c("max_wt_loss_pct", "wt_loss_slope_max", "wt_rmssd_rate", "wt_acf1")) %>% 
  mutate(variable = factor(variable, levels = c("max_wt_loss_pct", "wt_loss_slope_max", "wt_rmssd_rate", "wt_acf1"))) %>%
  ggplot(aes(x = wt_loss_high, y = rcvm, fill = variable)) +
  ## Amacero colors from the poisonfrogs R package
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
           fill = "#C6E74B", alpha = 0.3) +
  geom_hline(yintercept = 25, linetype = 'dashed',linewidth = 1.2, color = '#C6E74B') +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(labels = c(max_wt_loss_pct = "Max Loss %", 
                               wt_loss_slope_max = "Max Slope", 
                               wt_rmssd_rate = "RMSSD", 
                               wt_acf1 = "ACF1"),
                    values = c("max_wt_loss_pct" = "#14374A", 
                               "wt_loss_slope_max" = "#3981BE", 
                               "wt_rmssd_rate" = "#64B936", 
                               "wt_acf1" = "#7D2A1D")) +
  labs(x = "Weight Loss High", y = "Coefficient of Variation - Weight", fill = element_blank()) +
  #facet_grid(~HA) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom',
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = alpha("#C6E74B", 0.3), color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))


summarize_rcvm(data = fullMorbidity, 
               vars = c("temp_5", "temp_rmssd_rate", "temp_acf1")) %>%
  mutate(variable = factor(variable, levels = c("temp_5", "temp_rmssd_rate", "temp_acf1"))) %>%
  ggplot(aes(x = wt_loss_high, y = rcvm, fill = variable)) +
  ## Amacero colors from the poisonfrogs R package - Peru
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
           fill = "#C6E74B", alpha = 0.3) +
  geom_hline(yintercept = 25, linetype = 'dashed',linewidth = 1.2, color = '#C6E74B') +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(labels = c(temp_5 = "Temp Day 1-5", 
                               temp_rmssd_rate = "RMSSD", 
                               temp_acf1 = "ACF1"),
                    values = c("temp_5" = "#3981BE", 
                               "temp_rmssd_rate" = "#64B936", 
                               "temp_acf1" = "#7D2A1D")) +
  labs(x = "Weight Loss High", y = "Coefficient of Variation - Temp", fill = element_blank()) +
  #facet_grid(~HA) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom',
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = alpha("#C6E74B", 0.3), color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))

 
########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
