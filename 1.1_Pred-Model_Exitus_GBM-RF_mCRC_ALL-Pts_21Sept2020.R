library(WriteXLS)
library(gdata)
library(dplyr)
library(ggplot2)

#install.packages("doParallel",lib="/home/slahoz/R")
#library("doParallel",lib.loc="/home/slahoz/R")

#library("foreach")   # Package "foreach" is necessary to then load package "glmnet" --- Even though, package "foreach" is a DEPENDENCY of package "doParallel", so it is NOT necessary to upload it (actually, if we upload it parallely to "doParallel", then there're some functions that get masked)!!!
#library("glmnet")   # Package "glmnet" is necessary to run Radua's function "simple_imputation"

#install.packages(c("crossval", "randomForest"))
library(crossval)
library(randomForest)



###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### I) PREDICTOR VARIABLES (21st September 2020) ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

    ############
    ### (1) IMPORT Input dataframe:

  ### Input dataframe:
setwd("~/Desktop/Analysis/Metastatic-mCRC_JMaurel_2020/Clinical-Survival_Data")

dat.mv <- read.csv2("Survival-Multivariate-R_mCRC_200pts-QT_26May2020.csv",
                    header=TRUE,sep=";",stringsAsFactors=FALSE)
dim(dat.mv)
#[1] 200  35


  ### X co-variables (PREDICTORS):

# Sex, Age, Number-Organs_Affected, ECOG_PS, LDH, ALP (Alkalyne Phosph.),
# Leukos, PCR, MSI, Location_Right-Left, Mutated-Genes (TP53, SMAD4 & FBXW7, 
# plus/not BRAF, PIK3CA)


    ############
    ### (2) FIX Input dataframe & Clinical variables:

### In Mutation-Status columns, 1 means "Mutated" and 0 means "Wild-Type":
#dat.mv[,colnames(dat.mv) %in% c("KRAS","BRAF","FBXW7","NRAS","PIK3CA","PTEN","SMAD4","TP53")] <- ifelse(dat.mv[,colnames(dat.mv) %in% c("KRAS","BRAF","FBXW7","NRAS","PIK3CA","PTEN","SMAD4","TP53")]==1,"Mut","No_Mut")


### Create column "RAS" (which includes KRAS/NRAS):
dat.mv$RAS <- 0

for(i in 1:nrow(dat.mv))
{
  if(dat.mv$KRAS[i]==1 | dat.mv$NRAS[i]==1)
  { dat.mv$RAS[i] <- 1 }
}


### CLASS of column-variables:
sapply(dat.mv,class)

dat.mv$Sex <- as.character(dat.mv$Sex)
#dat.mv$Sex <- ifelse(dat.mv$Sex=="H","Man","Woman")
dat.mv$Sex <- ifelse(dat.mv$Sex=="H",1,0)

all.aug$Exitus <- factor(ifelse(all.aug$Exitus==0,"Alive","Exitus"),
                         levels=c("Alive","Exitus"))

colnames(all.aug)[18] <- "Age"
#colnames(all.aug) <- make.names(colnames(all.aug), unique=TRUE)


# We want this variable "Num_Lesions_Hep_1is1to3_2is4to9_3isMore9" to be CATEGORICAL:
dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9 <- as.character(dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9)
dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9 <- gsub("9","3",dat.mv$Num_Lesions_Hep_1is1to3_2is4to9_3isMore9)

dat.mv$Number_Organs <- as.character(dat.mv$Number_Organs)

dat.mv$Number_Organs_2 <- dat.mv$Number_Organs
dat.mv$Number_Organs_2 <- gsub(0,1,dat.mv$Number_Organs_2)
dat.mv$Number_Organs_2 <- ifelse(dat.mv$Number_Organs_2==1,0,1) # 0 means one-organ-affected, and 1 means More-than-1!

dat.mv$Number_Organs <- as.numeric(dat.mv$Number_Organs)
dat.mv$Number_Organs_2 <- as.numeric(dat.mv$Number_Organs_2)

#dat.mv$ECOG_PS_DxM1 <- as.character(dat.mv$ECOG_PS_DxM1)

dat.mv$Alk_Phosph_.116 <- ifelse(dat.mv$Alk_Phosph>=116,1,0)

dat.mv$LEUCOS_.11100 <- ifelse(dat.mv$LEUCOS>=11100,1,0)

dat.mv$PCR <- as.numeric(as.character(dat.mv$PCR))
dat.mv$PCR_.1 <- as.numeric(dat.mv$PCR_.1)
dat.mv$PCR_.1 <- ifelse(dat.mv$PCR>=1,1,0)

dat.mv$CEA <- as.numeric(as.character(dat.mv$CEA))
dat.mv$CEA_.5 <- ifelse(dat.mv$CEA>=5,1,0)

dat.mv$MSS0_MSI1 <- as.numeric(as.character(dat.mv$MSS0_MSI1))

dat.mv$GEMCAD <- as.character(dat.mv$GEMCAD)

dat.mv$Location_Right1_Left2_Transverse3 <- as.character(dat.mv$Location_Right1_Left2_Transverse3)
dat.mv$Location_Right1_Left2_Transverse3 <- gsub("3","1",dat.mv$Location_Right1_Left2_Transverse3)
colnames(dat.mv)[34] <- "Location_Right1_Left2"
dat.mv$Location_Right1_Left2 <- as.numeric(dat.mv$Location_Right1_Left2)


## FIX "Treatment" variable:
levels(factor(dat.mv$Treatment))
#[1] "CAPOX"        "FOLFIRI"      "FOLFIRI+BEV"  "FOLFIRI+CET" 
#[5] "FOLFIRI+PANI" "FOLFOX"       "FOLFOX+BEV"   "FOLFOX+CET"  
#[9] "FOLFOX+PANI"  "XELOX"        "XELOX+BEV"

for(i in 1:nrow(dat.mv))
{
  if(dat.mv$Treatment[i] %in% c("CAPOX","FOLFIRI","FOLFOX","XELOX"))
  { dat.mv$Treatment_Group[i] <- "Doublet" }
  
  else 
  { dat.mv$Treatment_Group[i] <- "Doublet_plus_Targeted-Agent" }
}

#
nrow(dat.mv[dat.mv$Treatment_Group=="Doublet",])
#[1] 139
nrow(dat.mv[dat.mv$Treatment_Group=="Doublet_plus_Targeted-Agent",])
#[1] 61



    ############
    ### (3) Generate FINAL-Dataframe for the ML-model:

### X co-variables (PREDICTORS):

# Sex, Age, Number-Organs_Affected, ECOG_PS, LDH, ALP (Alkalyne Phosph.),
# Leukos, PCR, MSI, Location_Right-Left, Mutated-Genes (TP53, SMAD4 & FBXW7, 
# plus/not BRAF, PIK3CA)


### COMBINED Model (Clínico-Mutacional):
all.pfs <- dat.mv[,c(2,11,15:18,37,21:23,25,27,31,34)] # To predict Progression
all.os <- dat.mv[,c(3,11,15:18,37,21:23,25,27,31,34)] # To predict Exitus (v1)
all.os <- dat.mv[,c(3,9:11,13,15:18,37,21:23,25,27,31,34)] # To predict Exitus (v2)
  
### GENETICO-MUTATIONAL Model:
#all.pfs <- dat.mv[,2:4]
#all.os <- dat.mv[,2:4]

### CLINICAL Model:
#all.pfs <- dat.mv[,c(1,5:11)]
#all.os <- dat.mv[,c(1,5:11)]

### GEMCAD-only Model:
#all.pfs <- dat.mv[,c(1,5:11)]
#all.os <- dat.mv[,c(1,5:11)]



###############################################
######### To ADD the Days_to_Recurrence and Days_to_Exitus within "all_var2" dataframe:

# (???)



###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
### ###### ###### II) PIPELINE --- PREDICTIVE MODEL (21st September 2020) ###### ###
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

https://topepo.github.io/caret/model-training-and-tuning.html

# Overall USED PACKAGES: "mice", "caret", "crossval", "pROC"

#all.aug <- all.pfs
all.aug <- all.os


####################################
####################################
### (2) [NO CAL!!!] IMPUTATION of Missing Values [package "mice"]:

library(mice)
# https://www.analyticsvidhya.com/blog/2016/03/tutorial-powerful-packages-imputing-missing-values/

## METHODS of Imputation by MICE package (function "mice"):
# i- PMM (Predictive Mean Matching)  – For NUMERIC Variables
# ii- logreg (Logistic Regression) – For BINARY Variables (with 2 levels)
# iii- polyreg (Bayesian polytomous regression) – For FACTOR Variables (>= 2 levels)
# iv- Proportional odds model (ordered, >= 2 levels)


##################
### 2.0) CHECK the Missing-Values in our data/Variables:

## mice package has a function known as md.pattern() -- It returns a TABULAR-Form (TABLE) of missing values 
## present in each variable in the data set:
md.pattern(all.aug)

# RESULT: Only "Number_Subclonal_CNAs" shows Missing NA-Values!


## Create a PLOT + TABLE which represents missing values:
#install.packages("VIM",dependencies=TRUE)
library(VIM)

aggr(all.aug, col=c('navyblue','yellow'),
     numbers=TRUE, sortVars=TRUE,
     labels=names(all.aug), cex.axis=.7,
     gap=3, ylab=c("Missing data","Pattern"))



##################
### 2.1) Imputation of NUMERICAL Variables: Age, CD8_quantitative_intratumor, Number_Subclonal_CNAs, CNA_Load

    ### Use imputation-METHOD "pmm" (only consider IMPUTATIONS/PREDICTIONS for "Number_Subclonal_CNAs"):
imputed_Data <- mice(all.aug, m=50, maxit = 500, method = 'pmm', seed = 500) # Possible METHODS of imputation***
# *** https://www.rdocumentation.org/packages/mice/versions/2.25/topics/mice
# Check METHODS "pmm" or "rf"!
summary(imputed_Data)


    ### CHECK imputed values in ECOG_PS_DxM1 -- Concretely, we check 
    ### the MEAN of ECOG_PS_DxM1 for the 1st patient/row that shows
    ### an NA-value for this variable ECOG_PS_DxM1:
dim(imputed_Data$imp$ECOG_PS_DxM1) # N=2 patients have been IMPUTED!
#[1] 2 50
imputed_Data$imp$ECOG_PS_DxM1  # Check which rows/PATIENTS were NA for ECOG_PS_DxM1!
mean(unlist(imputed_Data$imp$ECOG_PS_DxM1[1,])) # 1st patient
#[1] 0.4


## CHECK imputed values in LDH --
dim(imputed_Data$imp$LDH) # N=16 patients have been IMPUTED!
#[1] 16 50
#imputed_Data$imp$LDH  # Check which rows/PATIENTS were NA for LDH!
mean(unlist(imputed_Data$imp$LDH[1,]))
#[1] 365.52   # This 1st imputed-patient have been predicted the 40% of times as "0" (= wild-type), and the 60% of times as "1" (= mutated)!


## CHECK imputed values in Alk_Phosph --
dim(imputed_Data$imp$Alk_Phosph) # N=14 patients have been IMPUTED!
#[1] 14 50
#imputed_Data$imp$Alk_Phosph  # Check which rows/PATIENTS were NA for Alk_Phosph!
mean(unlist(imputed_Data$imp$Alk_Phosph[1,]))
#[1] 92.14   # This 1st imputed-patient have been predicted the 40% of times as "0" (= wild-type), and the 60% of times as "1" (= mutated)!


## CHECK imputed values in LEUCOS --
dim(imputed_Data$imp$LEUCOS) # N=11 patients have been IMPUTED!
#[1] 11 50
#imputed_Data$imp$LEUCOS  # Check which rows/PATIENTS were NA for LEUCOS!
mean(unlist(imputed_Data$imp$LEUCOS[1,]))
#[1] 6820.8   # This 1st imputed-patient have been predicted the 40% of times as "0" (= wild-type), and the 60% of times as "1" (= mutated)!


## CHECK imputed values in PCR --
dim(imputed_Data$imp$PCR) # N=53 patients have been IMPUTED!
#[1] 53 50
#imputed_Data$imp$PCR  # Check which rows/PATIENTS were NA for PCR!
mean(unlist(imputed_Data$imp$PCR[1,]))
#[1] 6.1116   # This 1st imputed-patient have been predicted the 40% of times as "0" (= wild-type), and the 60% of times as "1" (= mutated)!


## CHECK imputed values in MSS0_MSI1 --
dim(imputed_Data$imp$MSS0_MSI1) # N=1 patient has been IMPUTED!
#[1] 1 50
imputed_Data$imp$MSS0_MSI1  # Check which rows/PATIENTS were NA for MSS0_MSI1!
mean(unlist(imputed_Data$imp$MSS0_MSI1[1,])) # 1st patient
#[1] 0   # This 1st imputed-patient have been predicted the 40% of times as "0" (= wild-type), and the 60% of times as "1" (= mutated)!



    ### (A) Get COMPLETE-Data for ECOG_PS_DxM1/MSS0_MSI1 -- Obtain the MEAN value of the
    ### computed 50 pmm-matrixes/models:

## 1- ECOG_PS_DxM1 -- Get COMPLETE Data for NA-cases (n=2) of ECOG_PS_DxM1:
      # Dataframe for ECOG_PS_DxM1:
completeData.ecog <- data.frame(matrix(nrow=nrow(imputed_Data$imp$ECOG_PS_DxM1), ncol=1))
colnames(completeData.ecog)[1] <- "Mean_Imputed_Value"
rownames(completeData.ecog) <- rownames(imputed_Data$imp$ECOG_PS_DxM1)

      # Loop-ECOG_PS_DxM1:
for(i in 1:nrow(imputed_Data$imp$ECOG_PS_DxM1))
{ completeData.ecog[i,"Mean_Imputed_Value"] <- mean(unlist(imputed_Data$imp$ECOG_PS_DxM1[i,])) }

for(i in 1:nrow(completeData.ecog)){ 
  if(completeData.ecog$Mean_Imputed_Value[i] <0.5) 
    { completeData.ecog$Mean_Imputed_Value[i] <- 0 }
  else if(completeData.ecog$Mean_Imputed_Value[i] >=0.5 && completeData.ecog$Mean_Imputed_Value[i] <1.5) 
    { completeData.ecog$Mean_Imputed_Value[i] <- 1 }
  else if(completeData.ecog$Mean_Imputed_Value[i] >=1.5 && completeData.ecog$Mean_Imputed_Value[i] <2.5)
    { completeData.ecog$Mean_Imputed_Value[i] <- 2 }
  else if(completeData.ecog$Mean_Imputed_Value[i] >=2.5) 
    { completeData.ecog$Mean_Imputed_Value[i] <- 3 }
}

      # Add the MEAN imputed-ECOG_PS_DxM1 for each of n=6 CASES in the dataframe "all.aug":
all.aug[rownames(completeData.ecog),"ECOG_PS_DxM1"] <- completeData.ecog[,"Mean_Imputed_Value"]


## 2- MSS0_MSI1 -- Get COMPLETE Data for the NA-case (n=1) of MSS0_MSI1:
      # Dataframe for MSS0_MSI1:
completeData.msi <- data.frame(matrix(nrow=nrow(imputed_Data$imp$MSS0_MSI1), ncol=1))
colnames(completeData.msi)[1] <- "Mean_Imputed_Value"
rownames(completeData.msi) <- rownames(imputed_Data$imp$MSS0_MSI1)

      # Loop-MSS0_MSI1:
for(i in 1:nrow(imputed_Data$imp$MSS0_MSI1))
{ completeData.msi[i,"Mean_Imputed_Value"] <- mean(unlist(imputed_Data$imp$MSS0_MSI1[i,])) }

      # Add the MEAN imputed-MSS0_MSI1 for the n=1 CASE in the dataframe "all.aug":
all.aug[rownames(completeData.msi),"MSS0_MSI1"] <- completeData.msi[,"Mean_Imputed_Value"]


## Checkings:
all.aug[is.na(all.aug$ECOG_PS_DxM1),]
completeData.ecog
all.aug[rownames(completeData.ecog),]



    ### (B) Get COMPLETE-Data for LDH/Alk_Phosph/LEUCOS/PCR -- Obtain the MEAN value of the
    ### computed 50 pmm-matrixes/models:

## 1- LDH -- Get COMPLETE Data for NA-cases (n=16) of LDH cases:
        # Dataframe for LDH:
completeData.ldh <- data.frame(matrix(nrow=nrow(imputed_Data$imp$LDH), ncol=1))
colnames(completeData.ldh)[1] <- "Imputed_LDH"
rownames(completeData.ldh) <- rownames(imputed_Data$imp$LDH)

        # Loop-LDH:
for(i in 1:nrow(imputed_Data$imp$LDH))
{ completeData.ldh[i,"Imputed_LDH"] <- mean(unlist(imputed_Data$imp$LDH[i,])) }

        # Add the MEAN imputed-LDH for each of n=16 CASES in the dataframe "all.aug":
all.aug[rownames(completeData.ldh),"LDH"] <- completeData.ldh[,"Imputed_LDH"]


## 2- Alk_Phosph:
      # Dataframe for Alk_Phosph:
completeData.alkp <- data.frame(matrix(nrow=nrow(imputed_Data$imp$Alk_Phosph), ncol=1))
colnames(completeData.alkp)[1] <- "Imputed_Alk_Phosph"
rownames(completeData.alkp) <- rownames(imputed_Data$imp$Alk_Phosph)

      # Loop-Alk_Phosph:
for(i in 1:nrow(imputed_Data$imp$Alk_Phosph))
{ completeData.alkp[i,"Imputed_Alk_Phosph"] <- mean(unlist(imputed_Data$imp$Alk_Phosph[i,])) }

      # Add the MEAN imputed-Alk_Phosph for each of n=14 CASES in the dataframe "all.aug":
all.aug[rownames(completeData.alkp),"Alk_Phosph"] <- completeData.alkp[,"Imputed_Alk_Phosph"]


## 3- LEUCOS:
      # Dataframe for LEUCOS:
completeData.leucos <- data.frame(matrix(nrow=nrow(imputed_Data$imp$LEUCOS), ncol=1))
colnames(completeData.leucos)[1] <- "Imputed_LEUCOS"
rownames(completeData.leucos) <- rownames(imputed_Data$imp$LEUCOS)

      # Loop-LEUCOS:
for(i in 1:nrow(imputed_Data$imp$LEUCOS))
{ completeData.leucos[i,"Imputed_LEUCOS"] <- mean(unlist(imputed_Data$imp$LEUCOS[i,])) }

      # Add the MEAN imputed-LEUCOS for each of n=11 CASES in the dataframe "all.aug":
all.aug[rownames(completeData.leucos),"LEUCOS"] <- completeData.leucos[,"Imputed_LEUCOS"]


## 4- PCR:
      # Dataframe for PCR:
completeData.pcr <- data.frame(matrix(nrow=nrow(imputed_Data$imp$PCR), ncol=1))
colnames(completeData.pcr)[1] <- "Imputed_PCR"
rownames(completeData.pcr) <- rownames(imputed_Data$imp$PCR)

      # Loop-PCR:
for(i in 1:nrow(imputed_Data$imp$PCR))
{ completeData.pcr[i,"Imputed_PCR"] <- mean(unlist(imputed_Data$imp$PCR[i,])) }

      # Add the MEAN imputed-PCR for each of n=53 CASES in the dataframe "all.aug":
all.aug[rownames(completeData.pcr),"PCR"] <- completeData.pcr[,"Imputed_PCR"]


## Checkings:
all.aug[is.na(all.aug$PCR),]
completeData.pcr
all.aug[rownames(completeData.pcr),]


#
secure.imputed <- all.aug
secure <- all.aug

#write.table(secure,"~/Desktop/Secure-Checkpoint_Pred-Model_mCRC_After-Imputation_All-Variables_22Sept2020_v3.txt",
#             col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



####################################
####################################
### (3.0) Possible 3 MODELS:

# a- CLINICO-GENOMIC Model:
#all.aug <- secure[,c(1:10,12:20)] # Remove row 10 ("CDX2_intratumor") and 21 ("KRAS_Mut")!
#all.aug <- secure # Consider ALL variables!

# b- CLINICAL Model (or Clinico-Pathological):
#all.aug <- secure[,c(1:10,12:13)]   # Clinico-Pathological Model (X-variables)

# c- GENOMIC Model:
#all.aug <- secure[,c(1,14:20)] # Genomic Model (X-variables, without KRAS_Mut)



####################################
####################################
### (3.1) STANDARDIZATION/RE-CENTERING of variables 
### (We have added this standardization-step within the "train" function):

#nums <- unlist(lapply(all.aug, is.numeric))  

#st.all.aug <- data.frame(scale(all.aug[,nums]))
#st.all.aug <- cbind(data.frame(Status=all.aug[,"Status"]),st.all.aug)



####################################
####################################
### (3.2) [GOOD] CROSS-VALIDATION:

#all <- all.aug
#all <- st.all.aug


###### ###### ######
### For CROSS-VALIDATION:
n.folds = 4

    ##
folds = data.frame(Fold=1:200 %% n.folds)

#number_folds <- 0:(n.folds - 1)   # For Leave-One-Out crossvalidation
number_folds <- c(0,1,2,3)   # For k-fold=10 cross-validation


###### ###### ######
### MODELS:

    # (a) CLINICO-GENETIC [21st September]!
#all <- cbind(folds, all.aug)
#all <- cbind(folds, all.aug[,c(1,4,6:17)])
#all <- cbind(folds, st.all.aug)

    # (b) CLINICAL Model:
#all <- cbind(folds, all.aug[,c(1,8:17)])

    # (c) GENETIC Model (v1):
#all <- cbind(folds, all.aug[,1:7])

    # (d) GENETIC Model (v2):
all <- cbind(folds, all.aug[,c(1,4,6:7)])

    # (e) GEMCAD-variables Model:
gemcad <- cbind(all.aug,kkk1)
all <- cbind(folds,gemcad)
  

###### ###### ######
### FIX the Response-Variable (Y):
all.aug$Exitus <- factor(ifelse(all.aug$Exitus==1,"Exitus","Alive"),levels=c("Exitus","Alive"))


###### ###### ######
### FIT the ML-Classifier:
for(f in 1:length(number_folds))
{
  
  fold <- number_folds[f]
  
  print(paste0("Fold Number: ",fold))
  
  ###
  
  train <- all[all$Fold != fold,]   # All (k=10) TRAINING-datasets must have 82-83 individuals --- In this case, "train1" contains individuals/rows that do NOT have a 0 in column "Fold"
  train <- train[,-1]
  
  test <- all[all$Fold == fold,]   # All (k=10) TEST-datasets must have the remaining 10 individuals (not used in the paired training-datasets)
  test <- test[,-1]
  
  ###
  
  #test.pred.probs <- data.frame(row.names=c(1:nrow(test)))
  
  
  ####################################
  ####################################
  ### (4) FIT the CLASSIFIER MODEL (in this case, "GBM") [package "caret"]:
  
  ### (3.2) Internal CROSS-VALIDATION [package "caret"]:
  #library(caret)
  
  # Determine Training-Test:
  
  # Type of cross-validation (k-fold with 10 folds):
  #set.seed(321) # For the Clinico-Genomic Model and the Clinical Model!
  
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated 100 times
    repeats = 100,
    classProbs=TRUE)
  
  
  
  ################### GBM (Stochastic Gradient Boosting Machine):
  #gbmFit1 <- train(Exitus ~ ., data = train, preProcess=c("center","scale"), 
  #                 method = "gbm", # Initial method I tested: "gbm"
  #                 trControl = fitControl,
  #                 ## This last option is actually one
  #                 ## for gbm() that passes through
  #                 verbose = FALSE)
  
  
  
  ################### RandomForest:
  
  library(randomForest)
  library(mlbench)
  library(e1071)
  
  #https://rpubs.com/phamdinhkhanh/389752
  
  mtry <- sqrt(ncol(train))
  tunegrid <- expand.grid(.mtry=mtry)
  
  gbmFit1 <- train(Exitus ~ ., data = train, 
                   preProcess=c("center","scale"), 
                   method = "rf",
                   metric="Accuracy", 
                   tuneGrid = tunegrid, 
                   trControl = fitControl,
  ## This last option is actually one
  ## for gbm() that passes through
                   verbose = FALSE)
  
  
  
  ####################################
  ####################################
  ### (5) Extracting PREDICTIONS and PROBABILITIES (Numeric or Class/Categoric) [package "caret"]:
  
  #test.preds <- predict(gbmFit1, data=train, newdata = test)
  test.pred.probs <- predict(gbmFit1, data=train, newdata = test, type = "prob")
  rownames(test.pred.probs) <- rownames(test)
  
  #test.preds2 <- data.frame(class=ifelse(test.preds=="Recurrent","1","0"))
  #test.preds2$class <- as.factor(test.preds2$class)
  
  
  if(fold == 0)
  { prob.folds <- test.pred.probs }
  
  else if (fold != 0)
  { prob.folds <- rbind(prob.folds,test.pred.probs) }
  # 10 Folds
}

dim(prob.folds)
#[1] 200  2

#seccc <- prob.folds



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
#[1] 0.285


# CONFUSION MATRIX (False-Positives, False-Negative, True-Positives, True-Negatives):
cm <- crossval::confusionMatrix(pred.class2, real.class2, negative="0")
cm   
#FP TP TN FN 
#29 71 74 26


# PERFORMANCE PARAMETERS (Accuracy, Sensitivity, Specificity, PPV, NPV):


### (B) GOOD!!! CLINICO-GENETIC MODEL [21st September] --- RandomForest (adding KRAS, BRAF and PIK3CA):
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#0.7550000 0.7802198 0.7339450 0.7100000 0.8000000 2.2816784

### (B.2) GOOD!!! CLINICO-GENETIC MODEL [21st September] --- RandomForest (without KRAS, BRAF and PIK3CA):
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#0.7400000 0.7500000 0.7307692 0.7200000 0.7600000 2.0971411

### (B.3) GOOD!!! CLINICAL MODEL [21st September] --- RandomForest: kkk1
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#0.7500000 0.7604167 0.7403846 0.7300000 0.7700000 2.2029338 

### (B.4) GOOD!!! GENETIC MODEL [21st September] --- RandomForest (adding KRAS, BRAF and PIK3CA):
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
#0.5700000 0.5875000 0.5583333 0.4700000 0.6700000 0.5880407

### (B.5) GOOD!!! GENETIC MODEL [21st September] --- RandomForest (without KRAS, BRAF and PIK3CA):
diagnosticErrors(cm)
#acc      sens      spec       ppv       npv       lor 
# 0.4750000  0.4468085  0.4836601  0.2100000  0.7400000 -0.2789569



####################################
####################################
### (6.5) Get CONFIDENCE INTERVALS of Sensitivity, Specificity, PPV & NPV:
library(bdpv)
#https://rdrr.io/cran/bdpv/man/BDtest.html
#http://estadistica-dma.ulpgc.es/doctoMed/Tarea_6_-_Curvas_ROC.html # How to APPLY/USE BDtest function!

ci.mat <- as.matrix(data.frame(True_positive = cm[c(2,4)],
                               True_negative = cm[c(1,3)],
                               row.names=c("Test_positive","Test_negative")))

# In our case, the PREVALENCE is the (observed/real) PROPORTION of Recurrent-patients in our COHORT:
prop.table(table(all.aug$Exitus))["Exitus"]
#Exitus 
#0.5

# Apply BDtest function:
ci.values <- BDtest(ci.mat, pr=0.5, conf.level=0.90) # The prevalence is the REAL-proportion of Recurrent-patients!

# Outputs of BDtest function:
ci.values$INDAT
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
ci.values$PPVNPVDAT # Estimations of PPV & NPV, & their 95% CI


###### ###### ######
### (B) GOOD!!! CLINICO-GENETIC MODEL [21st September] --- RandomForest (adding KRAS, BRAF and PIK3CA):
ci.values$INDAT
#                 True positive True negative
#Test positive            71            29
#Test negative            20            80
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.7802198       0.7146283       0.6968790       0.8493339
#Specificity 0.7339450       0.6721143       0.6554168       0.8026560
ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.7695560       0.7195268       0.7041967       0.8240816
#PPV 0.7457121       0.7026369       0.6896872       0.7946339


###### ###### ######
### (B) GOOD!!! CLINICO-GENETIC MODEL [21st September] --- RandomForest (without KRAS, BRAF and PIK3CA):
#              True positive True negative
#Test positive            72            28
#Test negative            24            76
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.7500000       0.6844721       0.6667985       0.8212314
#Specificity 0.7307692       0.6671315       0.6499773       0.8014021
ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.7450980       0.6971139       0.6826174       0.7989022
#PPV 0.7358491       0.6908563       0.6773568       0.7870706


###### ###### ######
### (B) GOOD!!! CLINICAL MODEL [21st September] --- RandomForest:
ci.values$INDAT
#               True positive True negative
#Test positive            73            27
#Test negative            23            77
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.7604167       0.6955088       0.6779739       0.8303764
#Specificity 0.7403846       0.6772352       0.6601884       0.8099894
ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.7555192       0.7075785       0.6930266       0.8087997
#PPV 0.7454831       0.7005894       0.6870679       0.7962280


###### ###### ######
### (B) GOOD!!! GENETIC MODEL [21st September] --- RandomForest (adding KRAS, BRAF and PIK3CA):
ci.values$INDAT
#               True positive True negative
#Test positive            47            53
#Test negative            33            67
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.5875000       0.5096956       0.4894077       0.6805912
#Specificity 0.5583333       0.4956923       0.4791448       0.6353211
ci.values$PPVNPVDAT 
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.5751073       0.5256164       0.5114509       0.6363665
#PPV 0.5708502       0.5267819       0.5141822       0.6257197


###### ###### ######
### (B) GOOD!!! GENETIC MODEL [21st September] --- RandomForest (without KRAS, BRAF and PIK3CA):
ci.values$INDAT
#               True positive True negative
#Test positive            21            79
#Test negative            26            74
ci.values$SESPDAT # Estimations of Sensitivity & Specificity, & their 95% CI
#             Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#Sensitivity 0.4468085       0.3466438       0.3222928       0.5765589
#Specificity 0.4836601       0.4290327       0.4146252       0.5531759
ci.values$PPVNPVDAT
#     Estimate Lower 90% limit Lower 95% limit Upper 95% limit
#NPV 0.4664700       0.4173820       0.4037167       0.5303026
#PPV 0.4639041       0.4072019       0.3915017       0.5378609


####################################
####################################
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
roc.obj1 <- roc(real.class3 ~ pred.probs22) # CLINICO-GENETIC Model [21st September] -- Adding KRAS, BRAF and PIK3CA
roc.obj2 <- roc(real.class3 ~ pred.probs22) # CLINICO-GENETIC Model [21st September] -- Without KRAS, BRAF and PIK3CA
roc.obj_B <- roc(real.class3 ~ pred.probs22) # CLINICAL Model [21st September]
roc.obj_C1 <- roc(real.class3 ~ pred.probs22) # GENETIC Model [21st September] -- Adding KRAS, BRAF and PIK3CA
roc.obj_C2 <- roc(real.class3 ~ pred.probs22) # GENETIC Model [21st September] -- Without KRAS, BRAF and PIK3CA


coords <- coords(roc.obj, x="best", input="threshold") # OR "x="all""!
#  threshold specificity sensitivity
#1     0.609         0.9        0.64


auc(real.class3 ~ pred.probs22) # GBM - CLINICO-GENETIC Model [21st September]
#Area under the curve: 0.7766
auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September]
#Area under the curve: 0.8169

auc.aug1 <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September] --- Adding KRAS, BRAF and PIK3CA
#Area under the curve: 0.8169
auc.aug2 <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICO-GENETIC Model [21st September] --- Without KRAS, BRAF and PIK3CA
#Area under the curve: 0.8128
auc.aug_B <- auc(real.class3 ~ pred.probs22) # RandomForest - CLINICAL Model [21st September]
#Area under the curve: 0.8124
auc.aug_C1 <- auc(real.class3 ~ pred.probs22) # RandomForest - GENETIC Model [21st September] --- Adding KRAS, BRAF and PIK3CA
#Area under the curve: 0.5648 
auc.aug_C2 <- auc(real.class3 ~ pred.probs22) # RandomForest - GENETIC Model [21st September] --- Without KRAS, BRAF and PIK3CA
#Area under the curve: 0.5462


ci.auc(auc.aug1) # CLINICO-GENETIC Model --- Test DeLong for == 95% Confidence-Intervals (CIs) of AUC:
#95% CI: 0.7577-0.8761 (DeLong)
ci.auc(auc.aug2)
#95% CI: 0.7531-0.8725 (DeLong)
ci.auc(auc.aug_B) # CLINICAL Model [21st September]
#95% CI: 0.7522-0.8726 (DeLong)
ci.auc(auc.aug_C1) # GENETIC Model [21st September]
#95% CI: 0.4844-0.6451 (DeLong)
ci.auc(auc.aug_C2)
#95% CI: 0.4665-0.6259 (DeLong)


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


### BAD!!!
library(ROC)
ROC::ROC(pred.probs22,real.class3)



  