library(glmnet)
library(caret)


####################################
####################################
### (3.2) [GOOD] CROSS-VALIDATION:
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Predictive-Model_ML/Checkpoint_files")

all.aug <- read.table("Secure-Checkpoint_Pred-Model_mCRC_After-Imputation_All-Variables_21Sept2020_v2.txt",
                      header=TRUE,stringsAsFactors=FALSE)
#all.aug <- read.table("Secure-Checkpoint_Pred-Model_DEATH_mCRC_After-Imputation_All-Variables_22Apr2021.txt",
#                  header=TRUE,stringsAsFactors=FALSE)
#all <- all.aug
#all <- st.all.aug



###### ###### ######
### FIX the Response-Variable (Y):
all.aug$Exitus <- ifelse(all.aug$Exitus==1,"Exitus","Alive")#,levels=c("Exitus","Alive"))



###### ###### ######
### For CROSS-VALIDATION: kkk1
n.folds = 5

##
folds = data.frame(Fold=1:200 %% n.folds)

#number_folds <- 0:(n.folds - 1)   # For Leave-One-Out crossvalidation
number_folds <- c(0,1,2,3,4)   # For k-fold=10 cross-validation


###### ###### ######
### MODELS:

# (a) CLINICO-GENETIC [19th February]!
#all <- cbind(folds, all.aug)
#all <- cbind(folds, all.aug[,c(1,4,6:7)])
#all <- cbind(folds, st.all.aug)

# (b) CLINICAL Model:
#all <- cbind(folds, all.aug[,c(1,8:17)])
#all <- cbind(folds, all.aug[,c(1,5:ncol(all.aug))])

# (d) GENETIC Model (v1):
#all <- cbind(folds, all.aug[,c(1,4,6:7)])

# (c) GENETIC Model (v2):
all <- cbind(folds, all.aug[,1:7])
#all <- cbind(folds, all.aug[,1:4])

# (e) GEMCAD-variables Model:
#gemcad <- cbind(all.aug,kkk1)
#all <- cbind(folds,gemcad)



###### ###### ######
### FIT the ML-Classifier:
for(f in 1:length(number_folds))
{
  
  fold <- number_folds[f]
  
  print(paste0("Fold Number: ",fold))
  
  ###
  
  train <- all[all$Fold != fold,]   # All (k=10) TRAINING-datasets must have 82-83 individuals --- In this case, "train1" contains individuals/rows that do NOT have a 0 in column "Fold"
  #train <- train[,-1]
  
  test <- all[all$Fold == fold,]   # All (k=10) TEST-datasets must have the remaining 10 individuals (not used in the paired training-datasets)
  #test <- test[,-1]
  
  ###
  
  #test.pred.probs <- data.frame(row.names=c(1:nrow(test)))
  
  
  
  
  ####################################
  ####################################
  ### (7) FEATURE-SELECTION by CV LASSO-based Importance-analysis [package glmnet]:
  
  #https://www.machinelearningplus.com/machine-learning/feature-selection/#3lassoregression
  
  trainData <- train[,c(2:ncol(train))] # Perform the LASSO-based feature-selection on the Training-set!
  
  #
  x <- data.matrix(trainData[,-1],rownames.force=NA) # all X vars
  #y <- as.double(as.matrix(trainData[, 1])) # Only Class
  y <- as.double(as.matrix(ifelse(trainData[, 1]=='Exitus', 1, 0))) # Only Class
  
  
  # Fit the LASSO model on the TRAINING-dataset (Lasso: Alpha = 1):
  set.seed(982743)
  cv.lasso <- cv.glmnet(x, y, family='binomial', #type.measure='class',
                        nfolds=4, alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')
  
  # Results:
  #plot(cv.lasso)
  
  # plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
  cat('Min Lambda: ', cv.lasso$lambda.min, '\n 1Sd Lambda: ', cv.lasso$lambda.1se)
  df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
  df_coef2 <- data.frame(Variable=rownames(df_coef),Coef=df_coef)
  
  # See all CONTRIBUTING-VARIABLES (with Beta-coef > 0) + the Intercept:
  #df_coef[df_coef[, 1] != 0, ]
  beta.vars <- as.character(df_coef2$Variable[df_coef2$s1!=0 & df_coef2$Variable!="(Intercept)"])
  
  
  ####################################
  ####################################
  ### (8) FIT the CLASSIFIER MODEL (in this case, "GBM") [package "caret"]:
  
  ############
  ###### (8.1) Set seed for each Model:
  
  # CLINICO-GENOMIC Model:
  #set.seed(43879287) # 20 may 2022
  
  # CLINICAL Model:
  #set.seed(4398782) # 20 may 2022
  
  # GENOMIC Model:
  set.seed(43879287) # 20 may 2022
  #set.seed(573869)
  
  
  
  ############
  ###### (4) Internal CROSS-VALIDATION (k-fold with 10 folds):
  #library(caret)
  
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    repeats = 100, ## repeated 100 times
    classProbs=TRUE)
  
  
  
  ############
  ###### (5) Run the classifier GBM (Stochastic Gradient Boosting Machine) algorithm:
  
  ## VERSION 1 -- GOOD (library "caret"): kkk1
  gbmFit1 <- train(Exitus ~ ., data = train[,colnames(train) %in% c(beta.vars,"Exitus")], # Only include those variables present in "beta.vars" (pre-select by LASSO)!
                   preProcess=c("center","scale"),
                   method = "gbm", # Initial method I tested: "gbm"
                   trControl = fitControl,
                   #metric=c("Accuracy","RMSE"),
                   ## This last option is actually one
                   ## for gbm() that passes through
                   verbose = FALSE)
  
  
  ####################################
  ####################################
  ### (6) Extracting PREDICTIONS and PROBABILITIES (Numeric or Class/Categoric) [package "caret"]:
  
  #test.preds <- predict(gbmFit1, data=train, newdata = test)
  test.pred.probs <- stats::predict(gbmFit1, 
                                    data=train[,colnames(train) %in% c(beta.vars,"Exitus")],
                                    newdata = test,#[,colnames(test) %in% c(beta.vars,"Exitus")],
                                    type = "prob") # type="response" (if we want to get Numerical-PROBABILITIES!)
  rownames(test.pred.probs) <- rownames(test)
  
  #test.preds2 <- data.frame(class=ifelse(test.preds=="Recurrent","1","0"))
  #test.preds2$class <- as.factor(test.preds2$class)
  
  
  if(fold == 0) { 
    prob.folds <- test.pred.probs
    sel.vars <- data.frame(Variables=beta.vars,
                           Fold=rep(fold,length(beta.vars)))
  }
  
  else if (fold != 0) {
    prob.folds <- rbind(prob.folds,test.pred.probs) 
    sel.vars <- rbind(sel.vars,data.frame(Variables=beta.vars,
                                          Fold=rep(fold,length(beta.vars))))
  }
  # 10 Folds
}

dim(prob.folds)
#[1] 200  1
dim(sel.vars)
#[1] 52  2



### 19th May 2022:
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Tables_Paper/Supplementary_Tables")
WriteXLS(sel.vars,"Table-S6_Feature-Selection_Pred-Model_mCRC_May2022.xls",
         col.names=TRUE)



### Keep the probabilities from the 3 models:
clin.gen <- prob.folds # CLINICO-GENETIC Model
#clin <- prob.folds # CLINICAL Model
#gen <- prob.folds # GENETIC Model



####################################
####################################
### (6) EVALUATE the PREDICTIVE CAPACITY of the Trained-Model [package "crossval"]:

### EVALUATE the Accuracy/Sensitivity/Specificity of the Model:
library("crossval")


## 6.1- PREDICTED Classes (derived from our Statistical-Prediction) using 0.5 cut-off:
prob.folds$Pred_class <- colnames(prob.folds)[apply(prob.folds,1,which.max)]

pred.class <- data.frame(class=ifelse(prob.folds$Pred_class=="Exitus", "1", "0"))   # We can use the fitted prob to get) the prediction
rownames(pred.class) <- rownames(prob.folds)


## 6.2- REAL Classes (derived from Real-Observations):
real.class <- data.frame(class=all[order(match(rownames(all),rownames(pred.class))),"Exitus"])
rownames(real.class) <- rownames(pred.class)
real.class$class <- ifelse(real.class$class=="Exitus","1","0")


# COMPLETE MODEL (Clinico-Genomic):

## Training error rate (= fraction of missclassified cases):
pred.class2 <- factor(pred.class$class,levels=c("0","1"))
real.class2 <- factor(real.class$class,levels=c("0","1"))

mean(pred.class2 != real.class2)    # We compare the prediction with the real-observations
#[1] 0.255


# CONFUSION MATRIX (False-Positives, False-Negative, True-Positives, True-Negatives):
cm <- crossval::confusionMatrix(pred.class2, real.class2, negative="0")
cm   
#FP TP TN FN 
#33 67 78 22 


# PERFORMANCE PARAMETERS (Accuracy, Sensitivity, Specificity, PPV, NPV):


### CLINICO-GENETIC MODEL [19th May 2022] --- GBM:
diagnosticErrors(cm)
#     acc      sens      spec       ppv       npv       lor 
#0.7500000 0.7777778 0.7272727 0.7000000 0.8000000 


### CLINICAL MODEL:
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#0.7300000 0.7446809 0.7169811 0.7000000 0.7600000 


### GENETIC MODEL_v1:
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 


### GENETIC MODEL_v2:
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#



####################################
####################################
### (6.5) Get CONFIDENCE INTERVALS of Accuracy:
d.err <- data.frame(diagnosticErrors(cm))

## 95% CIs in Clinico-Genetic Model:
d.err["acc",] + 
  qnorm(c(0.025,0.5,0.975)) * sqrt( ((1-d.err["acc",]) * d.err["acc",]) / 100)
#[1] 0.6374848 0.7250000 0.8125152


## 95% CIs in Clinical Model:
d.err["acc",] + 
  qnorm(c(0.025,0.5,0.975)) * sqrt( ((1-d.err["acc",]) * d.err["acc",]) / 100)
#[1] 0.6762932 0.7600000 0.8437068


## 95% CIs in Genetic Model v1:
d.err["acc",] + 
  qnorm(c(0.025,0.5,0.975)) * sqrt( ((1-d.err["acc",]) * d.err["acc",]) / 100)
#


## 95% CIs in Genetic Model v1:
d.err["acc",] + 
  qnorm(c(0.025,0.5,0.975)) * sqrt( ((1-d.err["acc",]) * d.err["acc",]) / 100)
#



####################################
####################################
### (6.6) Get CONFIDENCE INTERVALS of Sensitivity, Specificity, PPV & NPV:
library(bdpv)
#https://rdrr.io/cran/bdpv/man/BDtest.html
#http://estadistica-dma.ulpgc.es/doctoMed/Tarea_6_-_Curvas_ROC.html # How to APPLY/USE BDtest function!

ci.mat <- as.matrix(data.frame(True_positive = cm[c(2,4)],
                               True_negative = cm[c(1,3)],
                               row.names=c("Test_positive","Test_negative")))

# In our case, the PREVALENCE is the (observed/real) PROPORTION of Recurrent-patients in our COHORT:
#prop.table(table(all.aug$Exitus))["Exitus"]
prop.table(table(all.aug$Exitus))
#0   1 
#0.5 0.5 

# Apply BDtest function:
ci.values <- BDtest(ci.mat, pr=0.5, conf.level=0.95) # The prevalence is the REAL-proportion of Recurrent-patients!

# Outputs of BDtest function:
ci.values$INDAT
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
ci.values$PPVNPVDAT # Estimations of PPV & NPV, & their 95% CI


###### ###### ######
### (A) CLINICO-GENETIC MODEL [19th May 2022] --- GBM:
ci.values$INDAT
#              True positive True negative
#Test positive            67            33
#Test negative            22            78
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.7802198       0.6968790         0.6811797         0.8602830
#Specificity 0.7339450       0.6554168         0.6407376         0.8140271
ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.7695560       0.7041967         0.6905146         0.8332844
#PPV 0.7457121       0.6896872         0.6782084         0.8031650


###### ###### ######
### (B) CLINICAL MODEL:
ci.values$INDAT
#               True positive True negative

ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit

ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit




###### ###### ######
### (C) GENETIC MODEL v1 [23rd April] --- LASSO (only TP53, SMAD4 & FBXW7):
ci.values$INDAT
#               True positive True negative

ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit

ci.values$PPVNPVDAT 
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit



###### ###### ######
### (C) GENETIC MODEL v2 [23rd April] --- LASSO (adding KRAS, BRAF and PIK3CA):
ci.values$INDAT
#               True positive True negative

ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit

ci.values$PPVNPVDAT 
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit



########################################################################
############### [[NO calculem l'AUC a 23rd April!!!]]
########################################################################

### (7) DISCRIMINATORY CAPACITY of the included-VARIABLES/PREDICTORS (AUC + ROC) [package "pROC"]:
library(pROC)

## INPUT-dataframes:
#
prob.folds$Real_class <- all[order(match(rownames(all),rownames(prob.folds))),
                             "Exitus"]
#
for(i in 1:nrow(prob.folds))
{ prob.folds$Pred_prob[i] <- prob.folds[i,colnames(prob.folds)[colnames(prob.folds)==prob.folds$Real_class[i]]] }
#pred.probs2 <- unlist(prob.folds$Pred_prob)


#
prob.folds$Pred_prob_exit <- prob.folds$Exitus
pred.probs22 <- prob.folds$Exitus


### CORRECT --- Score AUC (Area Under the Curve):
real.class3 <- as.numeric(as.character(real.class2))
pred.class3 <- as.numeric(as.character(pred.class2))

# CORRECT AUC!!!!!!!!!!!
library(pROC)

roc.obj1 <- roc(real.class3 ~ pred.probs22) # CLINICO-GENETIC Model_v1 [19th May 2022] -- Adding KRAS, BRAF and PIK3CA
roc.obj2 <- roc(real.class3 ~ pred.probs22) # CLINICO-GENETIC Model_v2
roc.obj_B <- roc(real.class3 ~ pred.probs22) # CLINICAL Model_v1
roc.obj_C1 <- roc(real.class3 ~ pred.probs22) # GENETIC Model_v2
roc.obj_C2 <- roc(real.class3 ~ pred.probs22) # GENETIC Model_v3


coords <- coords(roc.obj, x="best", input="threshold") # OR "x="all""!
#  threshold specificity sensitivity
#1     0.609         0.9        0.64


auc(real.class3 ~ pred.probs22) # GBM - CLINICO-GENETIC Model_v1 [19th May 2022]
#Area under the curve: 0.7808
auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September]
#Area under the curve:

auc.aug1 <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September] --- Adding KRAS, BRAF and PIK3CA
#Area under the curve:
auc.aug2 <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September] --- Without KRAS, BRAF and PIK3CA
#Area under the curve:
auc.aug_B <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICAL Model [21st September]
#Area under the curve:
auc.aug_C1 <- auc(real.class3 ~ pred.probs22) # RandomForest - GENETIC Model [21st September] --- Adding KRAS, BRAF and PIK3CA
#Area under the curve: 
auc.aug_C2 <- auc(real.class3 ~ pred.probs22) # RandomForest - GENETIC Model [21st September] --- Without KRAS, BRAF and PIK3CA
#Area under the curve:


ci.auc(auc.aug1) # CLINICO-GENETIC Model --- Test DeLong for == 95% Confidence-Intervals (CIs) of AUC:
#95% CI: 0.7154-0.8463 (DeLong)
ci.auc(auc.aug2)
#95% CI:
ci.auc(auc.aug_B) # CLINICAL Model [21st September]
#95% CI:
ci.auc(auc.aug_C1) # GENETIC Model [21st September]
#95% CI:
ci.auc(auc.aug_C2)
#95% CI:


# Test DeLong to comopare AUC derived from 2 Models:
roc.test(roc.obj1,roc.obj2, method=c("delong"))$p.value
#[1] 1
roc.test(roc.obj1,roc.obj_B, method=c("delong"))$p.value # CLINICO-GENOMIC vs. CLINICAL Models [21st September]
#[1] 0.5072972
roc.test(roc.obj2,roc.obj_B, method=c("delong"))$p.value # CLINICO-GENOMIC vs. CLINICAL Models [21st September]
#[1] 0.9394609
roc.test(roc.obj,roc.obj_C1, method=c("delong"))$p.value # CLINICO-GENOMIC vs. GENETIC Models [21st September]
#
roc.test(roc.obj_B,roc.obj_C2, method=c("delong"))$p.value # CLINICAL vs. GENETIC Models [21st September]
#
roc.test(roc.obj_C1,roc.obj_C2, method=c("delong"))$p.value
#[1] 0.7954523

### Plot 3 ROC curves (in 1 graph):
roc_A <- plot(roc.obj1, print.auc=TRUE,
              col="#1F56CC")
roc_A <- plot(roc.obj_B, print.auc = TRUE, 
              col = "#EDB234", print.auc.y = .35, add = TRUE)
roc_A <- plot(roc.obj_C, print.auc = TRUE, 
              col = "#CB1D18", print.auc.y = .2, add = TRUE)

##
roc_A <- plot(roc.obj1, print.auc=TRUE,
              col="#00139C", print.auc.y = .4)
roc_A <- plot(roc.obj2, print.auc=TRUE,
              col="#0321FE", print.auc.y = .35, add = TRUE)
roc_A <- plot(roc.obj_B, print.auc = TRUE, 
              col = "#A1CB1D", print.auc.y = .3, add = TRUE)
roc_A <- plot(roc.obj_C1, print.auc = TRUE, 
              col = "#E88406", print.auc.y = .25, add = TRUE)
roc_A <- plot(roc.obj_C2, print.auc=TRUE,
              col="#D00808", print.auc.y = .2, add = TRUE)


# ROC curve (smooth):
plot.roc(smooth(roc.obj),
         col="blue", print.auc=TRUE)


# 3 ROC curves using "ggroc" (from library "ggplot2"):
ggroc(list(Clinico_Genomic = roc.obj, Clinical = roc.obj_B, Genomic = roc.obj_C)) +
  theme_classic() + geom_abline(intercept = 1, slope = 1, col="#757575") +
  scale_color_manual(values=c("#1F56CC", "#EDB234","#CB1D18")) +
  ylab("Sensitivity") + xlab("Specificity")


## Receiver Operating Characteristics (ROC) Curve:
myRoc <- roc(real.class3 ~ pred.probs22, auc.polygon=TRUE, grid=TRUE, #smooth=TRUE,
             auc=TRUE, ci=TRUE, plot=TRUE)
plot(myRoc)


