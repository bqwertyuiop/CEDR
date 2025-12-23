###############load Packages
library(dplyr)
library(ggplot2)
library(h2o)
library(otrimle)
library(ica)
###############set all the parameters
DAE_G=3
SAE_G=2
PCA_G=3
ICA_G=4
methyhidden <- 10
mRNAhidden <- 10
miRNAhidden <- 90
###############load real datasets
setwd("~/CEDR/data/LGG")
load("LGG.Rdata")

######################
# run CEDR
######################

################# Dimensionality reduction and OTRIMLE clustering            
######DAE
h2o.init(max_mem_size = '100g')
methyinput_dropout_ratio <- 0.1
mRNAinput_dropout_ratio <- 0.1
miRNAinput_dropout_ratio <- 0.1
data1 <- methy_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = methyhidden,
  epochs = 10,
  l1 = 0.005,
  l2 = 0.1,
  activation = 'RectifierWithDropout',
  input_dropout_ratio = methyinput_dropout_ratio,
  hidden_dropout_ratios = 0.3,
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
head(ae1_codings,500,500)
DAEmethyfeatures <- attributes(ae1_codings)[["data"]]

data1 <- mRNA_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = mRNAhidden,
  epochs = 10,
  l1 = 0.005,
  l2 = 0.1,
  activation = 'RectifierWithDropout',
  input_dropout_ratio = methyinput_dropout_ratio,
  hidden_dropout_ratios = 0.3,
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
head(ae1_codings,500,500)
DAEmRNAfeatures <- attributes(ae1_codings)[["data"]]

data1 <- miRNA_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = miRNAhidden,
  epochs = 10,
  l1 = 0.005,
  l2 = 0.1,
  activation = 'RectifierWithDropout',
  input_dropout_ratio = methyinput_dropout_ratio,
  hidden_dropout_ratios = 0.3,
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
head(ae1_codings,500,500)
DAEmiRNAfeatures <- attributes(ae1_codings)[["data"]]

allDAEfeatures <- cbind(DAEmethyfeatures,DAEmRNAfeatures,DAEmiRNAfeatures)
##### OTRIMLE
x <- allDAEfeatures
set.seed(1)   
A2 <- try(otrimle(data=x, G=DAE_G,erc=1,npr.max=0.05),silent=TRUE)
if('try-error' %in% class(A2))          
{A2 <- 0 }
group_OTRIMLE_DAE <- head(A2[["cluster"]],500)

######SAE
data1 <- methy_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = methyhidden,
  l1 = 0.003,
  l2 = 0.1,
  epochs = 10,
  average_activation = 0.9,
  sparsity_beta = 0.05,
  activation = 'Rectifier',
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
ae1_codings1 <- head(ae1_codings,500,500)
SAEmethyfeatures <- ae1_codings1

data1 <- mRNA_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = mRNAhidden,
  l1 = 0.003,
  l2 = 0.1,
  epochs = 10,
  average_activation = 0.9,
  sparsity_beta = 0.05,
  activation = 'Rectifier',
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
head(ae1_codings,500,500)
SAEmRNAfeatures <- attributes(ae1_codings)[["data"]]

data1 <- miRNA_delcv
data2 <- as.data.frame(data1)
features <- as.h2o(data2)
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = miRNAhidden,
  l1 = 0.003,
  l2 = 0.1,
  epochs = 10,
  average_activation = 0.9,
  sparsity_beta = 0.05,
  activation = 'Rectifier',
  reproducible=TRUE,
  seed = 123
)
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
head(ae1_codings,500,500)
SAEmiRNAfeatures <- attributes(ae1_codings)[["data"]]

allSAEfeatures <- cbind(SAEmethyfeatures,SAEmRNAfeatures,SAEmiRNAfeatures)
####### OTRIMLE
x <- allSAEfeatures
set.seed(1)   
A2   <- try(otrimle(data=x, G=SAE_G,erc=1,npr.max=0.05),silent=TRUE)
if('try-error' %in% class(A2)) 
{A2 <- 0}
group_OTRIMLE_SAE <- head(A2[["cluster"]],500)

######PCA
methyPCAfeature <- 5
mRNAPCAfeature <- 5
miRNAPCAfeature <- 5

data <- methy_delcv
class(data) <- "numeric"
data1.pca <- prcomp(data, center=TRUE,scale. = TRUE) 
data2 <- data1.pca$x
methy_PCA  <- data2[,1:methyPCAfeature]

data <- mRNA_delcv
class(data) <- "numeric"
data1.pca <- prcomp(data, center=TRUE,scale. = TRUE)   
data2 <- data1.pca$x
mRNA_PCA  <- data2[,1:mRNAPCAfeature]

data <- miRNA_delcv
class(data) <- "numeric"
data1.pca <- prcomp(data, center=TRUE,scale. = TRUE) 
data2 <- data1.pca$x
miRNA_PCA  <- data2[,1:miRNAPCAfeature]

allPCAfeatures <- cbind(methy_PCA,mRNA_PCA,miRNA_PCA)
#####OTRIMLE
x <- allPCAfeatures
set.seed(1)   

A2    <- try(otrimle(data=x, G=PCA_G,erc=1,npr.max=0.05),silent=TRUE)
if('try-error' %in% class(A2))       
{
  A2 <- 0                             
}
group_OTRIMLE_PCA <- head(A2[["cluster"]],500)

######ICA
ICA_methy <- icaimax(methy_delcv, nc = 3)
ICA_mRNA <- icaimax(mRNA_delcv, nc = 3)
ICA_miRNA <- icaimax(miRNA_delcv, nc = 3)
ICAallfeatures <- cbind(ICA_methy[["S"]],ICA_mRNA[["S"]],ICA_miRNA[["S"]])
#####OTRIMLE
x <- ICAallfeatures
set.seed(1)   
A2    <- try(otrimle(data=x, G=ICA_G,erc=1,npr.max=0.05),silent=TRUE)
if('try-error' %in% class(A2))      
{
  A2 <- 0                               
}
group_OTRIMLE_ICA <- head(A2[["cluster"]],500)

ensemble2 <- group_OTRIMLE_PCA
Zero2 <- which(ensemble2==0)
ensemble3 <- group_OTRIMLE_DAE
Zero3 <- which(ensemble3==0)
ensemble4 <- group_OTRIMLE_SAE
Zero4 <- which(ensemble4==0)
ensemble5 <- group_OTRIMLE_ICA
Zero5 <- which(ensemble5==0)

C <- c(Zero2,Zero3,Zero4,Zero5)
C <- unique(C, fromLast = TRUE)
write.table(C,'C')
if(length(C) == 0){
  
  ensemble2C <- as.numeric(ensemble2)
  ensemble3C <- as.numeric(ensemble3)
  ensemble4C <- as.numeric(ensemble4)
  ensemble5C <- as.numeric(ensemble5)
  
  allensemble <- cbind(ensemble2C, ensemble3C, ensemble4C, ensemble5C)
  
} else {
  
  ensemble2C <- ensemble2[-C]
  ensemble3C <- ensemble3[-C]
  ensemble4C <- ensemble4[-C]
  ensemble5C <- ensemble5[-C]
  
  allensemble <- cbind(ensemble2C, ensemble3C, ensemble4C, ensemble5C)
}
write.table(allensemble,'allensemble')

####################################################################
#############################  Mixture Model for Clustering Ensembles
CEDRclustering <-function (Y, MAX = NULL, rep, SEED = 1) 
{
  LL <- NULL
  clusterresults <- NULL
  ensemblecluster <- NULL
  mm_k <- NULL
  BIC <- NULL
  numclus <- NULL
  I <- NULL
  N <- dim(Y)[1]
  H <- dim(Y)[2]
  set.seed(SEED)
  if (is.null(MAX)) {
    MAX <- max(Y[!is.na(Y)])
  }
  for (M in 2:MAX) {
    for (r in 1:rep) {
      assign(paste("CR", r, sep = ""), EM(Y, M))
      if (r == 1) {
        CR = get(paste("CR", r, sep = ""))
      }
      else {
        if (get(paste("CR", r, sep = ""))[["Log Likelihood"]] > 
            CR[["Log Likelihood"]]) {
          CR <- get(paste("CR", r, sep = ""))
        }
      }
    }
    i <- 1
    if (length(LL) > 0 && (CR[["Number of clusters"]] != 
                           M || CR[["Log Likelihood"]] < LL[length(LL)])) {
      repeat {
        CR <- EM(Y, M)
        i <- i + 1
        if ((CR[["Number of clusters"]] == M && CR[["Log Likelihood"]] > 
             LL[length(LL)]) || i > 3) {
          break
        }
      }
    }
    I <- c(I, i)
    if (length(LL) > 0 && (CR[["Number of clusters"]] != 
                           M || CR[["Log Likelihood"]] < LL[length(LL)])) {
      break
    }
    numclus <- c(numclus, CR[["Number of clusters"]])
    clusterresults <- c(clusterresults, list(CR[["Cluster results"]]))
    LL <- c(LL, CR[["Log Likelihood"]])
    BIC <- c(BIC, CR[["BIC"]])
    final_k_BIC <- numclus[which(BIC == min(BIC))]
  }
  
  BICresultsummary <- paste(
    "The final optimal k is", 
    numclus[which(BIC == min(BIC))], 
    "which has the lowest BIC of", 
    min(BIC)
  )
  
  BICcluster <- clusterresults[[which(BIC == min(BIC))]]
  
  result <- list(
    BIC_result_summary = BICresultsummary,
    BICcluster = BICcluster, 
    final_k_BIC = final_k_BIC
  )
  
  return(result)
}

allensembleCEDR <- allensemble
allensembleCEDR <- as.matrix(allensembleCEDR)
Y = allensembleCEDR
CEDRensemble <- try(CEDRclustering(Y, rep = 3,SEED = 123),silent=TRUE)

if ('try-error' %in% class(CEDRensemble)) {
  CEDRensemble <- NA
}
group_CEDR_BIC <- try(CEDRensemble[["BICcluster"]], silent = TRUE)
if ('try-error' %in% class(group_CEDR_BIC)) {
  group_CEDR_BIC <- NA
}
group <- group_CEDR_BIC