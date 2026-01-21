library(MMUPHin)

dataset <- read.csv("ML_Input_Dataset.csv", row.names = 1)
dataset <- t(dataset)
dataset_metadata <- read.csv("ML_Input_Dataset_Metadata.csv", row.names = 1)
dataset_metadata <- as.matrix(dataset_metadata)


fit_adjust_batch <- adjust_batch(feature_abd = dataset,
                                 batch = "BioProject",
                                 covariates = "SubGroup",
                                 data = dataset_metadata,
                                 control = list(verbose = FALSE))

write.csv(fit_adjust_batch[["feature_abd_adj"]], "ML_Input_Dataset_batch_corrected.csv" )





set.seed(101)
# list <- list('PRJEB12565', 'PRJNA299077', 'PRJNA316126','PRJNA439311','PRJNA377739', 'PRJNA543785', 'PRJNA667797', 
#              'PRJNA683885', "PRJNA795297", 'PRJNA934046', 'PRJNA944199', 'PRJNA1021628')
# 
# random_index <- sample(length(list), 1)
# selected_dataset <- list[[random_index]]
# selected_dataset

library(splitstackshape)

#Stratified random sampling (80-20%)
split <- read.csv("ML_Input_Dataset_batch_corrected.csv", header = TRUE)

set.seed(101)
split_80 <- stratified(split, c("BioProject", "SubGroup"), 0.8, select = NULL, replace = FALSE)
write.csv(split_80, file = "split_80_genus_level_all_stable_exacerbation_batch_corrected.csv")

split_20 <- read.csv("split_20_Runs_Only_all_stable_rxacerbation.csv", header = TRUE)
split_20_samples <- merge(split, split_20, by.x = "Run")
write.csv(split_20_samples, file = "split_20_genus_level_all_stable_exacerbation_batch_corrected.csv")






##Machine_Learning
# Input_Data
input_train <- read.csv("split_80_genus_level_all_stable_exacerbation_batch_corrected.csv", header = TRUE)
Train <- input_train [, -1:-2]
Train$SubGroup <- as.factor(Train$SubGroup)


Blind_input <- read.csv("split_20_genus_level_all_stable_exacerbation_batch_corrected.csv", header = TRUE)
Blind <- Blind_input [,  -1:-2]
Blind$SubGroup <- as.factor(Blind$SubGroup)


library(caret)
library(mltools)
library(data.table)

##SVM_L Model
set.seed(78945)
seeds_svm_l <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_svm_l[[i]] <- rep(78945, 1)
seeds_svm_l <- lapply(seeds_svm_l, as.vector)

fitControl_SVM_L <- trainControl(
  method = "CV",
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary
)

set.seed(78945)
svm_L_model <- train(
  SubGroup ~ ., 
  data = Train, 
  method = "svmLinear", 
  trControl = fitControl_SVM_L, 
  verbose = TRUE, 
  importance = TRUE
)

# Training results
svm_L_model_results <- svm_L_model$results
svm_L_model_pred <- svm_L_model$pred
svm_L_model_cf <- confusionMatrix(data = svm_L_model_pred$pred, reference = svm_L_model_pred$obs, mode = "everything")
svm_L_model_MC_summary <- multiClassSummary(svm_L_model_pred, lev = levels(svm_L_model_pred$pred))

# Compute MCC (training)
svm_L_model_pred_dt <- as.data.table(svm_L_model_pred)
svm_L_model_MCC <- mcc(preds = svm_L_model_pred_dt$pred, actuals = svm_L_model_pred_dt$obs)

# Print training metrics
svm_L_model_MC_summary
svm_L_model_cf
print(paste("Training MCC:", svm_L_model_MCC))

# Test predictions
prediction_svm_L_model <- predict.train(svm_L_model, Blind, probability = TRUE)
pred_svm_L_model_prob <- predict.train(svm_L_model, Blind, type = "prob")
pred_svm_L_model_prob$pred <- prediction_svm_L_model
pred_svm_L_model_prob$obs <- Blind$SubGroup

# Test metrics
svm_L_test_MC <- multiClassSummary(pred_svm_L_model_prob, lev = levels(pred_svm_L_model_prob$pred))
svm_L_test_cf <- confusionMatrix(data = prediction_svm_L_model, reference = Blind$SubGroup, mode = "everything")

# Compute MCC (test)
pred_svm_L_model_prob_dt <- as.data.table(pred_svm_L_model_prob)
svm_L_test_MCC <- mcc(preds = pred_svm_L_model_prob_dt$pred, actuals = pred_svm_L_model_prob_dt$obs)

# Print test metrics
svm_L_test_MC
svm_L_test_cf
print(paste("Test MCC:", svm_L_test_MCC))

# Save model
saveRDS(svm_L_model, file = "svm_L_model.rds")

#SVM-RBF MODEL
set.seed(78945)
seeds_rbf <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_rbf[[i]] <- rep(78945, 200)
seeds_rbf <- lapply(seeds_rbf, as.vector)

fitControl_SVM_R <- trainControl(
  method = "CV",
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_rbf
)

svm_grid_1 <- expand.grid(sigma = seq(0, 1, length = 10), C = seq(0, 10, length = 20))

set.seed(78945)
svm_R_model <- train(
  SubGroup ~ ., 
  data = Train,
  method = "svmRadial",
  trControl = fitControl_SVM_R,
  verbose = TRUE,
  tuneGrid = svm_grid_1,
  importance = TRUE
)

write.csv(svm_R_model$results, file = "SVM_R_parameters_tunelength_20.csv")

# svm_rbf model on specific parameter
set.seed(78945)
seeds_rbf_2 <- vector(mode = "list", length = 6) 
for(i in 1:6) seeds_rbf_2[[i]] <- rep(78945, 1)
seeds_rbf_2 <- lapply(seeds_rbf_2, as.vector)

fitControl_SVM_R_3 <- trainControl(
  method = "CV",
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_rbf_2
)

svm_grid_3 <- expand.grid(sigma = 1, C = 10)

set.seed(78945)
svm_R_model_1 <- train(
  SubGroup ~ ., 
  data = Train, 
  method = "svmRadial", 
  trControl = fitControl_SVM_R_3, 
  verbose = TRUE, 
  tuneGrid = svm_grid_3, 
  importance = TRUE
)

# TRAIN results
svm_R_model_results <- svm_R_model_1$results
svm_R_model_pred <- svm_R_model_1$pred
svm_R_model_cf <- confusionMatrix(data = svm_R_model_pred$pred, reference = svm_R_model_pred$obs, mode = "everything")
svm_R_model_MC_summary <- multiClassSummary(svm_R_model_pred, lev = levels(svm_R_model_pred$pred))

# Compute MCC on training predictions
svm_R_model_pred_dt <- as.data.table(svm_R_model_pred)
svm_R_model_MCC <-mcc(preds = svm_R_model_pred_dt$pred, actuals = svm_R_model_pred_dt$obs)

# TEST results
prediction_svm_R_model <- predict.train(svm_R_model_1, Blind, probability = TRUE)
pred_svm_R_model_prob <- predict.train(svm_R_model_1, Blind, type = "prob")
pred_svm_R_model_prob$pred <- prediction_svm_R_model
pred_svm_R_model_prob$obs <- Blind$SubGroup

svm_R_test_MC <- multiClassSummary(pred_svm_R_model_prob, lev = levels(pred_svm_R_model_prob$pred))
svm_R_test_cf <- confusionMatrix(data = prediction_svm_R_model, reference = Blind$SubGroup, mode = "everything")

# Compute MCC on test predictions
pred_svm_R_model_prob_dt <- as.data.table(pred_svm_R_model_prob)
svm_R_test_MCC <- mcc(preds = pred_svm_R_model_prob_dt$pred, actuals = pred_svm_R_model_prob_dt$obs)

# Print summaries
svm_R_model_MC_summary
print(paste("Training MCC:", svm_R_model_MCC))
svm_R_model_cf

svm_R_test_MC
print(paste("Test MCC:", svm_R_test_MCC))
svm_R_test_cf
# Save model
saveRDS(svm_R_model_1, file = "svm_R_model.rds")

##Random Forest
set.seed(78945)
seeds_rf <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_rf[[i]] <- rep(78945, 100)
seeds_rf <- lapply(seeds_rf, as.vector)

grid_RF <- expand.grid(mtry = sample(1:ncol(Train), 100))

fitControl_RF <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  seeds = seeds_rf
)

set.seed(78945)
RF_model <- train(
  SubGroup ~ .,
  data = Train,
  method = "rf",
  trControl = fitControl_RF,
  verbose = TRUE,
  tuneGrid = grid_RF,
  importance = TRUE
)

write.csv(RF_model$results, file = "RF_Model_Parameters_Genus.csv")


#RF model on specific parameters
set.seed(78945)
seeds_rf_2 <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_rf_2[[i]] <- rep(78945, 1)
seeds_rf_2 <- lapply(seeds_rf_2, as.vector)

grid_RF_2 <- expand.grid(mtry = 40)

fitControl_RF_2 <- trainControl(
  method = "CV",
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_rf_2
)

set.seed(78945)
RF_model_1 <- train(
  SubGroup ~ .,
  data = Train,
  method = "rf",
  trControl = fitControl_RF_2,
  verbose = TRUE,
  tuneGrid = grid_RF_2,
  importance = TRUE
)

# Training metrics
RF_pred <- RF_model_1$pred
RF_cf <- confusionMatrix(data = RF_pred$pred, reference = RF_pred$obs, mode = "everything")
RF_MC_summary <- multiClassSummary(RF_pred, lev = levels(RF_pred$pred))

# Compute MCC (training)
RF_pred_dt <- as.data.table(RF_pred)
RF_train_MCC <- mcc(preds = RF_pred_dt$pred, actuals = RF_pred_dt$obs)

# Print training results
RF_MC_summary
RF_cf
print(paste("Training MCC:", RF_train_MCC))

# Variable importance
RF_var <- varImp(RF_model_1)
RF_var_imp <- RF_var$importance

# Test predictions
prediction_RF <- predict.train(RF_model_1, Blind, probability = TRUE)
pred_RF_prob <- predict.train(RF_model_1, Blind, type = "prob")
pred_RF_prob$pred <- prediction_RF
pred_RF_prob$obs <- Blind$SubGroup

# Test metrics
RF_test_MC <- multiClassSummary(pred_RF_prob, lev = levels(pred_RF_prob$pred))
RF_test_cf <- confusionMatrix(data = prediction_RF, reference = Blind$SubGroup, mode = "everything")

# Compute MCC (test)
pred_RF_prob_dt <- as.data.table(pred_RF_prob)
RF_test_MCC <- mcc(preds = pred_RF_prob_dt$pred, actuals = pred_RF_prob_dt$obs)

# Print test results
RF_test_MC
RF_test_cf
print(paste("Test MCC:", RF_test_MCC))

# Save model
saveRDS(RF_model_1, file = "RF_model.rds")

##Naive bayes
set.seed(78945)
seeds_nb <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_nb[[i]] <- rep(78945, 18)
seeds_nb <- lapply(seeds_nb, as.vector)

fitControl_nb <- trainControl(
  method = "CV",
  p = 0.8,
  number = 5,
  classProbs = TRUE,
  savePredictions = "all",
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  seeds = seeds_nb
)

grid <- expand.grid(
  laplace = c(0, 0.5, 1.0),
  usekernel = c(TRUE, FALSE),
  adjust = c(0, 0.5, 1.0)
)

set.seed(78945)
nb_model <- train(
  SubGroup ~ ., 
  data = Train,
  method = "naive_bayes",
  trControl = fitControl_nb,
  verbose = TRUE,
  tuneGrid = grid,
  importance = TRUE
)

write.csv(nb_model$results, file = "nb_model_parameters.csv")

#NB Model with Fixed Parameter #
set.seed(78945)
seeds_nb_1 <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_nb_1[[i]] <- rep(78945, 1)
seeds_nb_1 <- lapply(seeds_nb_1, as.vector)

fitControl_nb_1 <- trainControl(
  method = "CV",
  p = 0.8,
  number = 5,
  classProbs = TRUE,
  savePredictions = "all",
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  seeds = seeds_nb_1
)

set.seed(78945)
nb_model_1 <- train(
  SubGroup ~ ., 
  data = Train,
  method = "naive_bayes",
  trControl = fitControl_nb_1,
  tuneGrid = expand.grid(laplace = 0, usekernel = FALSE, adjust = 0),
  verbose = TRUE,
  importance = TRUE
)

#TRAINING METRICS
nb_pred <- nb_model_1$pred
nb_cf <- confusionMatrix(data = nb_pred$pred, reference = nb_pred$obs, mode = "everything")
nb_MC_summary <- multiClassSummary(nb_pred, lev = levels(nb_pred$pred))

# Compute MCC (training)
nb_pred_dt <- as.data.table(nb_pred)
nb_train_MCC <- mcc(preds = nb_pred_dt$pred, actuals = nb_pred_dt$obs)

# Print training results
nb_MC_summary
nb_cf
print(paste("Training MCC:", nb_train_MCC))

# Variable importance
nb_var <- varImp(nb_model_1)
nb_var_imp <- nb_var$importance

#TEST METRICS
prediction_nb <- predict.train(nb_model_1, Blind, probability = TRUE)
pred_nb_prob <- predict.train(nb_model_1, Blind, type = "prob")
pred_nb_prob$pred <- prediction_nb
pred_nb_prob$obs <- Blind$SubGroup

nb_test_MC <- multiClassSummary(pred_nb_prob, lev = levels(pred_nb_prob$pred))
nb_test_cf <- confusionMatrix(data = prediction_nb, reference = Blind$SubGroup, mode = "everything")

# Compute MCC (test)
pred_nb_prob_dt <- as.data.table(pred_nb_prob)
nb_test_MCC <- mcc(preds = pred_nb_prob_dt$pred, actuals = pred_nb_prob_dt$obs)

# Print test results
nb_test_MC
nb_test_cf
print(paste("Test MCC:", nb_test_MCC))

#SAVE MODEL
saveRDS(nb_model_1, file = "nb_model.rds")


##glm
set.seed(78945)
seeds_lr <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_lr[[i]] <- rep(78945, 1)
seeds_lr <- lapply(seeds_lr, as.vector)

fitControl_glm <- trainControl(
  method = "CV",
  savePredictions = "all",
  p = 0.8,
  verboseIter = TRUE,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_lr
)

set.seed(78945)
lr_model <- train(
  SubGroup ~ ., 
  data = Train,
  method = "glm",
  trControl = fitControl_glm
)

#TRAINING METRICS
lr_pred <- lr_model$pred
lr_cf <- confusionMatrix(data = lr_pred$pred, reference = lr_pred$obs, mode = "everything")
lr_MC_summary <- multiClassSummary(lr_pred, lev = levels(lr_pred$pred))

# Compute MCC (training)
lr_pred_dt <- as.data.table(lr_pred)
lr_train_MCC <- mcc(preds = lr_pred_dt$pred, actuals = lr_pred_dt$obs)

# Print training results
lr_MC_summary
lr_cf
print(paste("Training MCC:", lr_train_MCC))

# Variable importance
lr_var <- varImp(lr_model)
lr_var_imp <- lr_var$importance

#TEST METRICS
prediction_lr <- predict.train(lr_model, Blind, probability = TRUE)
pred_lr_prob <- predict.train(lr_model, Blind, type = "prob")
pred_lr_prob$pred <- prediction_lr
pred_lr_prob$obs <- Blind$SubGroup

lr_test_MC <- multiClassSummary(pred_lr_prob, lev = levels(pred_lr_prob$pred))
lr_test_cf <- confusionMatrix(data = prediction_lr, reference = Blind$SubGroup, mode = "everything")

# Compute MCC (test)
pred_lr_prob_dt <- as.data.table(pred_lr_prob)
lr_test_MCC <- mcc(preds = pred_lr_prob_dt$pred, actuals = pred_lr_prob_dt$obs)

# Print test results
lr_test_MC
lr_test_cf
print(paste("Test MCC:", lr_test_MCC))

#SAVE MODEL
saveRDS(lr_model, file = "glm_model.rds")



##mlpML
#Hyperparameter Tuning
set.seed(78945)
seeds_mlpml <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_mlpml[[i]] <- rep(78945, 180)
seeds_mlpml <- lapply(seeds_mlpml, as.vector)

fitControl_mlpML <- trainControl(
  method = "CV",
  savePredictions = "all",
  verboseIter = TRUE,
  p = 0.8,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_mlpml
)

set.seed(78945)
mlpML_model_all <- train(
  x = Train[, colnames(Train) != "SubGroup"],
  y = Train$SubGroup,
  method = "mlpML",
  trControl = fitControl_mlpML,
  verbose = TRUE,
  tuneGrid = expand.grid(
    layer1 = c(10, 20, 30, 40, 50),
    layer2 = c(0, 10, 20, 30, 40, 50),
    layer3 = c(0, 10, 20, 30, 40, 50)
  ),
  importance = TRUE
)

write.csv(mlpML_model_all$results, file = "mlpML_model_parameters.csv")

#MLPML Model (Fixed Params)
set.seed(78945)
seeds_mlpml_1 <- vector(mode = "list", length = 6)
for (i in 1:6) seeds_mlpml_1[[i]] <- rep(78945, 1)
seeds_mlpml_1 <- lapply(seeds_mlpml_1, as.vector)

fitControl_mlpML_1 <- trainControl(
  method = "CV",
  savePredictions = "all",
  verboseIter = TRUE,
  p = 0.8,
  classProbs = TRUE,
  number = 5,
  summaryFunction = multiClassSummary,
  seeds = seeds_mlpml_1
)

set.seed(78945)
mlpML_model_1 <- train(
  x = Train[, colnames(Train) != "SubGroup"],
  y = Train$SubGroup,
  method = "mlpML",
  trControl = fitControl_mlpML_1,
  verbose = TRUE,
  tuneGrid = expand.grid(layer1 = 50, layer2 = 0, layer3 = 0),
  importance = TRUE
)

#TRAINING METRICS
mlpML_pred <- mlpML_model_1$pred
mlpML_cf <- confusionMatrix(data = mlpML_pred$pred, reference = mlpML_pred$obs, mode = "everything")
mlpML_MC_summary <- multiClassSummary(mlpML_pred, lev = levels(mlpML_pred$pred))

# Compute MCC (training)
mlpML_pred_dt <- as.data.table(mlpML_pred)
mlpML_train_MCC <- mcc(preds = mlpML_pred_dt$pred, actuals = mlpML_pred_dt$obs)

# Print training results
mlpML_MC_summary
mlpML_cf
print(paste("Training MCC:", mlpML_train_MCC))

#TEST METRICS
prediction_mlpML <- predict.train(mlpML_model_1, Blind, probability = TRUE)
pred_mlpML_prob <- predict.train(mlpML_model_1, Blind, type = "prob")
pred_mlpML_prob$pred <- prediction_mlpML
pred_mlpML_prob$obs <- Blind$SubGroup

mlpML_test_MC <- multiClassSummary(pred_mlpML_prob, lev = levels(pred_mlpML_prob$pred))
mlpML_test_cf <- confusionMatrix(data = prediction_mlpML, reference = Blind$SubGroup, mode = "everything")

# Compute MCC (test)
pred_mlpML_prob_dt <- as.data.table(pred_mlpML_prob)
mlpML_test_MCC <- mcc(preds = pred_mlpML_prob_dt$pred, actuals = pred_mlpML_prob_dt$obs)

# Print test results
mlpML_test_MC
mlpML_test_cf
print(paste("Test MCC:", mlpML_test_MCC))

#SAVE MODEL
saveRDS(mlpML_model_1, file = "mlpML_model.rds")

##Variable Importance


