#This Rscript allows to classify COPD stable and exacerbation stages using machine learning approaches
library(phyloseq)
library(file2meco)
library(microeco)

ps <- readRDS("PRJNA377739_Amplicon-16S_ps_object.rds")


# Create phyloseq object to meco object
meco_object <- phyloseq2meco(ps)
meco_object$tidy_dataset()

# Filter pollution
suppressMessages(meco_object$filter_pollution(taxa = c("mitochondria", "chloroplast")))
meco_object$tidy_dataset()

# Filter based on abundance and detection
suppressMessages(meco_object$filter_taxa(rel_abund = 0.0001, freq = 0.05))
meco_object$tidy_dataset()

#Machine learning
meco_object$cal_abund()
t2 <- trans_classifier$new(dataset = meco_object, y.response = "SubGroup", x.predictors = "all")

as <- as.data.frame(t2$data_feature)
as$Run <- rownames(as)

sam <- as.data.frame(meco_object$sample_table)
sam$Run <- rownames(sam)

feature_table <- merge(as, sam, by.x = "Run")

write.csv(feature_table, file = "PRJNA683885_feature_table_all.csv")


a <- read.csv("PRJEB12565_ML_Input_1.csv")
b <- read.csv("PRJNA299077_ML_Input_1.csv")
c <- read.csv("PRJNA316126_ML_Input_1.csv")
d <- read.csv("PRJNA377739_ML_Input_1.csv")
e <- read.csv("PRJNA439311_ML_Input_1.csv")
f <- read.csv("PRJNA543785_ML_Input_1.csv")
g <- read.csv("PRJNA667797_ML_Input_1.csv")
h <- read.csv("PRJNA795297_ML_Input_1.csv")
i <- read.csv("PRJNA934046_ML_Input_1.csv")
j <- read.csv("PRJNA944199_ML_Input_1.csv")
k <- read.csv("PRJNA1021628_ML_Input_1.csv")
l <- read.csv("PRJNA683885_ML_Input_1.csv")

library(dplyr)

merged_bp <- full_join(a, b, by = "Taxa")
merged_bp_1 <- full_join(merged_bp, c, by = "Taxa")
merged_bp_2 <- full_join(merged_bp_1, d, by = "Taxa")
merged_bp_3 <- full_join(merged_bp_2, e, by = "Taxa")
merged_bp_4 <- full_join(merged_bp_3, f, by = "Taxa")
merged_bp_5 <- full_join(merged_bp_4, g, by = "Taxa")
merged_bp_6 <- full_join(merged_bp_5, h, by = "Taxa")
merged_bp_7 <- full_join(merged_bp_6, i, by = "Taxa")
merged_bp_8 <- full_join(merged_bp_7, j, by = "Taxa")
merged_bp_9 <- full_join(merged_bp_8, k, by = "Taxa")
merged_bp_10 <- full_join(merged_bp_9, l, by = "Taxa")



merged_bp_10 = merged_bp_10 %>% replace(is.na(.), 0)

write.csv(merged_bp_10, file = "ML_Input_Merged.csv")

ml_input <- read.csv("ML_Input_Merged.csv")

#add the metadata here
metadata <- read.csv("Metadata_all_stable_exacerbation_samples_MDPD.csv")

v <- merge(ml_input, metadata, by.x = "Run")

write.csv(v, file = "Final_ml_input_stable_exacerbation_Include_Frequent_Infrequent.csv")

library(splitstackshape)

#Stratified random sampling (80-20%)
split <- read.csv("Final_ml_input_stable_exacerbation_all_stable_exacerbation.csv", header = TRUE)

set.seed(101)
split_80 <- stratified(split, c("BioProject", "SubGroup"), 0.8, select = NULL, replace = FALSE)
write.csv(split_80, file = "split_80_genus_level_all_stable_exacerbation.csv")

split_20 <- read.csv("split_20_Runs_Only_all_stable_rxacerbation.csv", header = TRUE)
split_20_samples <- merge(split, split_20, by.x = "Run")
write.csv(split_20_samples, file = "split_20_genus_level_all_stable_exacerbation.csv")


# Input_Data
input_train <- read.csv("split_80_genus_level_only_frequent.csv", header = TRUE)
Train <- input_train [, -1:-2]
Train$SubGroup <- as.factor(Train$SubGroup)


Blind_input <- read.csv("split_20_genus_level_only_frequent.csv", header = TRUE)
Blind <- Blind_input [,  -1:-2]
Blind$SubGroup <- as.factor(Blind$SubGroup)

library(caret)


#SVM-Linear model
set.seed(78945)
seeds_svm_l <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_svm_l[[i]] <- rep(11010, 1)
seeds_svm_l <- lapply(seeds_svm_l, as.vector)
fitControl_SVM_L <- trainControl(method = "CV", savePredictions = "all", p = 0.7, verboseIter = TRUE,  classProbs = TRUE, number = 5, summaryFunction = multiClassSummary)
set.seed(78945)
svm_L_model <- train(SubGroup ~ ., data = Train, method = "svmLinear", trControl = fitControl_SVM_L, verbose = TRUE, importance = TRUE)

svm_L_model_results <- svm_L_model$results
svm_L_model_pred <- svm_L_model$pred
svm_L_model_cf <- confusionMatrix(data = svm_L_model_pred$pred, reference = svm_L_model_pred$obs, mode = "everything")
svm_L_model_MC_summary <- multiClassSummary(svm_L_model_pred, lev = levels(svm_L_model_pred$pred))
svm_L_model_MC_summary
svm_L_model_cf
prediction_svm_L_model <- predict.train(svm_L_model, Blind, probability = TRUE)
pred_svm_L_model_prob <- predict.train(svm_L_model, Blind, type = "prob")
pred_svm_L_model_prob$pred <- prediction_svm_L_model
pred_svm_L_model_prob$obs <- Blind$SubGroup
svm_L_test_MC <- multiClassSummary(pred_svm_L_model_prob, lev = levels(pred_svm_L_model_prob$pred))
svm_L_test_cf <- confusionMatrix(data = prediction_svm_L_model, reference = Blind$SubGroup, mode = "everything")
svm_L_test_MC
svm_L_test_cf

saveRDS(svm_L_model, file = "svm_L_model.rds")


#SVM-RBF MODEL
set.seed(11010)
seeds_rbf <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_rbf[[i]] <- rep(11010, 200)
seeds_rbf <- lapply(seeds_rbf, as.vector)
fitControl_SVM_R <- trainControl(method = "CV", savePredictions = "all", p = 0.8, verboseIter = TRUE,  classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_rbf)
svm_grid_1 <- expand.grid(sigma = seq(0, 1, length = 10), C = seq(0, 10, length = 20))
set.seed(11010)
svm_R_model <- train(SubGroup ~ ., data = Train, method = "svmRadial", trControl = fitControl_SVM_R, verbose = TRUE, tuneGrid = svm_grid_1, importance = TRUE)
write.csv(svm_R_model$results, file = "SVM_R_parameters_tunelength_20.csv")

#svm_rbf model on specific parameter
set.seed(11010)
seeds_rbf_2 <- vector(mode = "list", length = 6) 
for(i in 1:6) seeds_rbf_2[[i]] <- rep(11010, 1)
seeds_rbf_2 <- lapply(seeds_rbf_2, as.vector)
fitControl_SVM_R_3 <- trainControl(method = "CV", savePredictions = "all", p = 0.8, verboseIter = TRUE,  classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_rbf_2)
svm_grid_3 <- expand.grid(sigma = 1, C = 9.47)
set.seed(11010)
svm_R_model_1 <- train(SubGroup ~ ., data = Train, method = "svmRadial", trControl = fitControl_SVM_R_3, verbose = TRUE, tuneGrid = svm_grid_3, importance = TRUE)


svm_R_model_results <-svm_R_model_1$results
svm_R_model_pred <- svm_R_model_1$pred
svm_R_model_cf <- confusionMatrix(data = svm_R_model_pred$pred, reference = svm_R_model_pred$obs, mode = "everything")
svm_R_model_MC_summary <- multiClassSummary(svm_R_model_pred, lev = levels(svm_R_model_pred$pred))
svm_R_model_MC_summary
svm_R_model_cf
svm_R_model_1_var <- varImp(svm_R_model_1)
svm_R_model_1_var_imp <- svm_R_model_1_var$importance
prediction_svm_R_model <- predict.train(svm_R_model_1, Blind, probability = TRUE)
pred_svm_R_model_prob <- predict.train(svm_R_model_1, Blind, type = "prob")
pred_svm_R_model_prob$pred <- prediction_svm_R_model
pred_svm_R_model_prob$obs <- Blind$SubGroup
svm_R_test_MC <- multiClassSummary(pred_svm_R_model_prob, lev = levels(pred_svm_R_model_prob$pred))
svm_R_test_cf <- confusionMatrix(data = prediction_svm_R_model, reference = Blind$SubGroup, mode = "everything")
svm_R_test_MC
svm_R_test_cf

saveRDS(svm_R_model_1, file = "svm_R_model.rds")


#RANDOM FOREST MODEL
set.seed(11010)
seeds_rf <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_rf[[i]] <- rep(11010, 100) 
seeds_rf <- lapply(seeds_rf, as.vector)
grid_RF <- expand.grid(mtry = sample(1:ncol(Train), 100))
fitControl_RF <- trainControl(method = 'cv', number = 5, savePredictions = "all", p = 0.8, verboseIter = TRUE,  classProbs = TRUE, summaryFunction = multiClassSummary, seeds = seeds_rf)
set.seed(11010)
RF_model <- train(SubGroup ~ ., data = Train, method = "rf", trControl = fitControl_RF, verbose = TRUE, tuneGrid = grid_RF, importance = TRUE)
write.csv(RF_model$results, file = "RF_Model_Parameters_Genus.csv")

#RF model on specific parameter
set.seed(11010)
seeds_rf_2 <- vector(mode = "list", length = 6) 
for(i in 1:6) seeds_rf_2[[i]] <- rep(11010, 1)
seeds_rf_2 <- lapply(seeds_rf_2, as.vector)
grid_RF_2 <- expand.grid(mtry = 260)
fitControl_RF_2 <- trainControl(method = "CV", savePredictions = "all", p = 0.8, verboseIter = TRUE,  classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_rf_2)
set.seed(11010)
RF_model_1 <- train(SubGroup ~ ., data = Train, method = "rf", trControl = fitControl_RF_2, verbose = TRUE, tuneGrid = grid_RF_2, importance = TRUE)

RF_pred <- RF_model_1$pred
RF_cf <- confusionMatrix(data = RF_pred$pred, reference = RF_pred$obs, mode = "everything")
RF_MC_summary <- multiClassSummary(RF_pred, lev = levels(RF_pred$pred))
RF_results <- RF_model_1$results
RF_MC_summary
RF_cf
RF_var <- varImp(RF_model_1)
RF_var_imp <- RF_var$importance
RF_var <- varImp(RF_model_mtry_146)
RF_var_imp <- RF_var$importance
prediction_RF <- predict.train(RF_model_1, Blind, probability = TRUE)
pred_RF_prob <- predict.train(RF_model_1, Blind, type = "prob")
pred_RF_prob$pred <- prediction_RF
pred_RF_prob$obs <- Blind$SubGroup
RF_test_MC <- multiClassSummary(pred_RF_prob, lev = levels(pred_RF_prob$pred))
RF_test_cf <- confusionMatrix(data = prediction_RF, reference = Blind$SubGroup, mode = "everything")
RF_test_MC
RF_test_cf

saveRDS(RF_model_1, file = "RF_model.rds")


#Naive bayes Model
set.seed(11010)
seeds_nb <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_nb[[i]] <- rep(11010, 18)
seeds_nb <- lapply(seeds_nb, as.vector)
fitControl_nb <- trainControl(method = "CV", p = 0.8, number = 5, classProbs = TRUE, savePredictions = "all", summaryFunction = multiClassSummary, verboseIter = TRUE, seeds = seeds_nb)
grid <- expand.grid(laplace=c(0,0.5,1.0), usekernel = c(TRUE, FALSE), adjust=c(0,0.5,1.0))
set.seed(11010)
nb_model <- train(SubGroup ~ ., data = Train, method = "naive_bayes", trControl = fitControl_nb, verbose = TRUE, tuneGrid = grid, importance = TRUE)
write.csv(nb_model$results, file = "nb_model_parametrs.csv")

#NB model at specific parameter
set.seed(11010)
seeds_nb_1 <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_nb_1[[i]] <- rep(11010, 1)
seeds_nb_1 <- lapply(seeds_nb_1, as.vector)
fitControl_nb_1 <- trainControl(method = "CV", p = 0.8, number = 5, classProbs = TRUE, savePredictions = "all", summaryFunction = multiClassSummary, verboseIter = TRUE, seeds = seeds_nb_1)
set.seed(11010)
nb_model_1 <- train(SubGroup ~ ., data = Train, method = "naive_bayes", trControl = fitControl_nb_1, tuneGrid = expand.grid(laplace=0, usekernel = FALSE, adjust=0), verbose = TRUE, importance = TRUE)

nb_pred <- nb_model_1$pred
nb_cf <- confusionMatrix(data = nb_pred$pred, reference = nb_pred$obs, mode = "everything")
nb_MC_summary <- multiClassSummary(nb_pred, lev = levels(nb_pred$pred))
nb_results <- nb_model_1$results
nb_MC_summary
nb_cf
nb_var <- varImp(nb_model_1)
nb_var_imp <- nb_var$importance
nb_var <- varImp(nb_model_1)
nb_var_imp <- nb_var$importance
prediction_nb <- predict.train(nb_model_1, Blind, probability = TRUE)
pred_nb_prob <- predict.train(nb_model_1, Blind, type = "prob")
pred_nb_prob$pred <- prediction_nb
pred_nb_prob$obs <- Blind$SubGroup
nb_test_MC <- multiClassSummary(pred_nb_prob, lev = levels(pred_nb_prob$pred))
nb_test_cf <- confusionMatrix(data = prediction_nb, reference = Blind$SubGroup, mode = "everything")
nb_test_MC
nb_test_cf

saveRDS(nb_model_1, file = "nb_model.rds")


#glm logistic regression Model
set.seed(11010)
seeds_lr <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_lr[[i]] <- rep(11010, 1)
seeds_lr <- lapply(seeds_lr, as.vector)
fitControl_glm <- trainControl(method = "CV", savePredictions = "all", p = 0.8, verboseIter = TRUE,  classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_lr)
set.seed(11010)
lr_model <- train(SubGroup ~ ., data = Train, method = "glm", trControl = fitControl_glm)
lr_pred <- lr_model$pred
lr_cf <- confusionMatrix(data = lr_pred$pred, reference = lr_pred$obs, mode = "everything")
lr_twoC_summary <- multiClassSummary(lr_pred, lev = levels(lr_pred$pred))
lr_results <- lr_model$results
lr_twoC_summary
lr_cf
lr_var <- varImp(lr_model)
lr_var_imp <- lr_var$importance
lr_var <- varImp(lr_model)
lr_var_imp <- lr_var$importance
prediction_lr <- predict.train(lr_model, Blind, probability = TRUE)
pred_lr_prob <- predict.train(lr_model, Blind, type = "prob")
pred_lr_prob$pred <- prediction_lr
pred_lr_prob$obs <- Blind$SubGroup
lr_test_twoC <- multiClassSummary(pred_lr_prob, lev = levels(pred_lr_prob$pred))
lr_test_cf <- confusionMatrix(data = prediction_lr, reference = Blind$SubGroup, mode = "everything")
lr_test_twoC
lr_test_cf

saveRDS(lr_model, file = "glm_model.rds")

#mlpml Model
set.seed(11010)
seeds_mlpml <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_mlpml[[i]] <- rep(11010, 180)
seeds_mlpml <- lapply(seeds_mlpml, as.vector)

fitControl_mlpML <- trainControl(method = "CV", savePredictions = "all", verboseIter = TRUE, p = 0.8, classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_mlpml)
set.seed(11010)
mlpML_model_all <- train(x = Train[ , colnames(Train) != "SubGroup"],
                     y = Train$SubGroup, method = "mlpML",
                     trControl = fitControl_mlpML, verbose = TRUE,
                     tuneGrid = expand.grid(layer1=c(10, 20, 30, 40, 50), layer2=c(0, 10, 20, 30, 40, 50), layer3 = c(0, 10, 20, 30, 40, 50)), importance = TRUE)


#mlpml model at specific parameter
set.seed(11010)
seeds_mlpml_1 <- vector(mode = "list", length = 6)
for(i in 1:6) seeds_mlpml_1[[i]] <- rep(11010, 1) #no of parameters (sigma, 5, c 5, then 11010, 25)
seeds_mlpml_1 <- lapply(seeds_mlpml_1, as.vector)

fitControl_mlpML_1 <- trainControl(method = "CV", savePredictions = "all", verboseIter = TRUE, p = 0.8, classProbs = TRUE, number = 5, summaryFunction = multiClassSummary, seeds = seeds_mlpml_1)
set.seed(11010)
mlpML_model_1 <- train(x = Train[ , colnames(Train) != "SubGroup"], 
                       y = Train$SubGroup, method = "mlpML", 
                       trControl = fitControl_mlpML_1, verbose = TRUE,
                       tuneGrid = expand.grid(layer1=20, layer2=0, layer3 = 20), importance = TRUE)


mlpML_pred <- mlpML_model_1$pred
mlpML_cf <- confusionMatrix(data = mlpML_pred$pred, reference = mlpML_pred$obs, mode = "everything")
mlpML_MC_summary <- multiClassSummary(mlpML_pred, lev = levels(mlpML_pred$pred))
mlpML_results <- mlpML_model_1$results
mlpML_MC_summary
mlpML_cf
prediction_mlpML <- predict.train(mlpML_model_1, Blind, probability = TRUE)
pred_mlpML_prob <- predict.train(mlpML_model_1, Blind, type = "prob")
pred_mlpML_prob$pred <- prediction_mlpML
pred_mlpML_prob$obs <- Blind$SubGroup
mlpML_test_MC <- multiClassSummary(pred_mlpML_prob, lev = levels(pred_mlpML_prob$pred))
mlpML_test_cf <- confusionMatrix(data = prediction_mlpML, reference = Blind$SubGroup, mode = "everything")
mlpML_test_MC
mlpML_test_cf

saveRDS(mlpML_model_1, file = "mlpml_model.rds")
