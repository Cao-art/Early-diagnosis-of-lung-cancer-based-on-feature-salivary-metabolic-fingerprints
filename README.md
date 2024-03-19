# Early-diagnosis-of-lung-cancer-based-on-feature-salivary-metabolic-fingerprints
#Data normalization

        library('vegan')
        
#Data normalization, with samples as rows and features as columns

        data =read.csv("test.csv",row.names = 1,check.names = F)
        x<-as.matrix(data[,-1])
        norm<- decostand(x, method="range",MARGIN=2)
        
###Correlation analysis

    require(psych)
    require(ggcor)
    cor <- corr.test(t(norm),method = "spearman")####spearman
    
####RSD

     QC = read.csv("QC.csv",row.names = 1,check.names = F)#读入数据,列名为特征，行名为样本
    RSD = function(d){
    rsd <- c()
    for (i in 2:ncol(d)){
    rsd[i-1]= sd(d[,i])/mean(d[,i])
    }
    rsd = as.data.frame(rsd)
    row.names(rsd) = colnames(d[,-1])
    return(rsd)
    }
    QCrsd = RSD(QC)
   
#### Normality test

    QC = read.csv("QC.csv",row.names = 1,check.names = F)
 re = c()
 for(i in 2:ncol(QC)){
  a = shapiro.test(as.numeric(QC[,i]))[["p.value"]]
  re = cbind(re,a)
 }
 re=data.frame(t(re))
 
### PCA

  rm(list =ls())
  data = read.csv("test.csv",check.names = F)
  data$Group<-as.factor(data$Group)
  norm1_pca<-prcomp(data[,-1])
  norm1_pcs<-data.frame(norm1_pca$x)
  norm1_pcs<-data.frame(norm1_pca$x,Group=data$Group)
 
#### Extracting the variance contribution rate of principal components

 summ1 <- summary(norm1_pca)
 xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
 ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

### PLS-DA

 library(ropls)
 plsda<- opls(x=data[,-1],y=data[,"Group"],permI=100,orthoI = 0)
 var_summary <- varianceSummary(plsda)
 
# Extracting the explained variance of T1 and T2 components

 x_lab <-paste0("t1(",round(plsda@modelDF[1, "R2X"] * 100,2),"%)")
 y_lab <-paste0("t1(",round(plsda@modelDF[2, "R2X"] * 100,2),"%)")
 
######Analysis of Variance (ANOVA)
##p-value

 library(plyr)
 library(reshape2)
 library(rstatix)
 wilcoxtest<- ddply(melt(data),"variable",
                   function(x){
                     w<- wilcox.test(value~Group,data=x,paired = F)##paired = c(T,F)
                     with(w,data.frame(statistic,p.value))
                   })
                   
##fdr

  FDR = adjust_pvalue(wilcoxtest,method = "fdr")
 
####FC

   FC = function(d){
   FC = c()
   for(i in 2:ncol(d)){
    FC[i-1] = mean(d[d$Group =="A",i])/mean(d[d$Group =="B",i])
    re = cbind(colnames(d[,-1]),FC)
    }
    }
   FC = as.data.frame(FC(data))  
   result = cbind(FDR,FC)[,-5]
   result$FC = as.numeric(result$FC)
   for(i in 1:nrow(result)){
   result[i,"log2FC"] = log2(result$FC[i])
  }
 
###metabolic functional analysis

    library(tidyverse)
    library(magrittr)
    library(clusterProfiler)
    library(reshape2)
    keggannotation <- read.csv("KEGG pathway enrichment.csv",check.names = F)##download from KEGG
    x <- clusterProfiler::enricher(gene = kegg_a$`kegg ID`,##metabolites KEGG ID
                               TERM2GENE =keggannotation,minGSSize = 1)
    re = as.data.frame(x@result) 
 
####machine learning

   rm(list =ls())
   data = read.csv("test.csv",check.names = F)
   library(caret)
   library(readxl)
   library(foreach)
   library(doParallel)
 
##feature selection---rfe

   cl <- makeCluster(detectCores()-1)
   registerDoParallel(cl)
   time1 <- Sys.time()#
   set.seed(0830)
   rfs = rfe(data[,-1], data[,1],sizes = seq(1,646,1),
          rfeControl = rfeControl(functions =rfFuncs,#caretFuncs
                                  saveDetails = T,
                                  method ="cv")
          #,method = "glm"
   )
   time2 <- Sys.time()
   print(time2-time1)
   stopImplicitCluster()
   stopCluster(cl)##
   dim(norm) #
   plot(rfs, type = c('g','o'))
 
###

   ctrl = trainControl(method = "cv",
                    number = 10,#fold or iterations
                    
                    search = "grid",###how the tuning parameter grid is determined
                    savePredictions = "final", 
                    returnResamp = "final",
                    summaryFunction =multiClassSummary,
                    classProbs =T,allowParallel = T,
                    verboseIter = T)
   cl <- makeCluster(detectCores()-1)
   registerDoParallel(cl)
   time1 <- Sys.time()#
   set.seed(1135)
   mods<- train(x=data[,-1],y=data[,1],method ='ranger',
             metric = "F1" ,#Selecting model evaluation metrics
             tuneLength =5,
             #tuneGrid = tune,
             trControl = ctrl
   )
   time2 <- Sys.time()
   print(time2-time1)
   stopImplicitCluster()
   stopCluster(cl)##
   saveRDS(mod1, "mod.rds")
   perfomance= function(m){
    Accuracy=mean(m[["resample"]][["Accuracy"]])
    Sensitivity=mean(m[["resample"]][["Sensitivity"]])
    Specificity=mean(m[["resample"]][["Specificity"]])
    Kappa= mean(m[["resample"]][["Kappa"]])
    F1=mean(m[["resample"]][["F1"]])
    Precision=mean(m[["resample"]][["Precision"]])
    Pos_Pred_Value=mean(m[["resample"]][["Pos_Pred_Value"]])
    Neg_Pred_Value=mean(m[["resample"]][["Neg_Pred_Value"]])
    re = data.frame(Accuracy,Sensitivity,Specificity,Pos_Pred_Value,Neg_Pred_Value,Kappa,F1,Precision)
    return(re)
    }
    reptest=perfomance(mods)###Outputting the results of model cross-validation

##model prediction

   prob<- predict(mods,newdata=test[,-1],type = "prob")
   pred<- predict(mods,newdata=test[,-1])
   confusionMatrix(data =prob$pred, reference = data$Group,positive = "A",mode = "everything")
 
####AUC\accuracy\sensitivity\specificity confidence interval

   library(boot)
   library(pROC)
   library(caret)
   roc_object <- roc(data$Group, prob$A,ci = T, levels=c("B","A"))
   cise<-ci.se(roc_object,conf.level=0.95,specificities=seq(0, 1, 0.01))
   prob$Group = data$Group
 
# Calculate the number of true positives, false positives, true negatives, and false negatives

   TP <- sum(prob$pred == "A" & prob$Group == "A")  # true positive
   FP <- sum(prob$pred == "A" & prob$Group == "B")  # false positive
   TN <- sum(prob$pred == "B" & prob$Group == "B")  # true negative
   FN <- sum(prob$pred == "B" & prob$Group == "A")  # false negative
 
# Calculate sensitivity, specificity, accuracy

   sensitivity <- TP / (TP + FN)
   specificity <- TN / (TN + FP)
   accuracy <- (TP + TN) / length(prob$Group)
 
# output results

   print(paste("Sensitivity:", sensitivity))
   print(paste("Specificity:", specificity))
   print(paste("Accuracy:", accuracy))
 
# Calculate the confidence interval for sensitivity using binom.test

   sensitivity_ci <- binom.test(x = TP, n = TP + FN, conf.level = 0.95)$conf.int
 
# Calculate the confidence interval for specificity using binom.test

   specificity_ci <- binom.test(x = TN, n = TN + FP, conf.level = 0.95)$conf.int
 
# Calculate the confidence interval for accuracy using binom.test.

   accuracy_ci <- binom.test(x = TP + TN, n = length(prob$Group), conf.level = 0.95)$conf.int

# output results

   print(paste("Sensitivity CI:", sensitivity_ci))
   print(paste("Specificity CI:", specificity_ci))
   print(paste("Accuracy CI:", accuracy_ci))

####ROC

   library(pROC)
   ROC =roc(prob$Group,prob$A,levels=c("B","A"))

##Chi-square test##

   a = c(79,13)
   b = c(66,4)
   c = c(65,4)
   s=chisq.test(cbind(a,cbind(b,c)),correct = F)
   s
   s$expected##expected frequency
   s$observed
 
###fisher.test

   fisher.test(cbind(a,cbind(b,c)))



