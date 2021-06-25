#######################################################
# SAMPLE CODE FOR GENERATING SIMULATED DATA FOR
# DYNAMIC SURVIVAL PREDICTION
#######################################################
####################################################
## ----------------------------------------------------------------------------
##  Code to accompany chapter 4 of the thesis by Kamaryn Tanner:
##    "Dynamic Prediction, Mediation and Communication for Survival
##      Outcomes, with Applications to Cystic Fibrosis"
##  
##  Contact for code: kamaryn.tanner1@lshtm.ac.uk
## ----------------------------------------------------------------------------
##  Date        :       24-June-2021
##  R Version   :       4.0.2                                                           
##  
##  Required files:     simParams.R -- contains parameters used in data generation           
##                      simSupportFxns.R -- contains helper functions called in this file
##  
##  R packages needed:  simjm, JM, plyr, dplyr, zoo, SuperLearner
##                      plus packages for algorithms used in SuperLearner
## ------------------------------------------------------------------
#####################################################################

## ---- Load R libraries ---- ##

# Information on simjm is available here: https://github.com/sambrilleman/simjm
library(simjm)
library(JM)
library(plyr)
library(dplyr)
library(SuperLearner)
library(ranger)  #for random forest algorithm
library(xgboost) #for boosting algorithm
library(zoo)

## ---- Load other files ---- ##

## Remember to set your working directory to the location of this file
source("simParams.R") #Contains parameters used to generate data,
                      # in the simjm call as well as for Super Learner and
                      # landmarking

source("simSupportFxns.R")  #Contains helper functions called in this file

## ---- Generate training, test data ---- ##

## set the seed
Seed <- 1212

## data is simulated from a joint model using simjm 
dataSim <- simjm(n = Ntrain*(R+1), M = 1, fixed_trajectory = "linear",
                     random_trajectory = "linear", assoc = "etavalue",
                     basehaz = c("weibull"), betaLong_intercept = L_intercept, betaLong_binary = L_binary,
                     betaLong_continuous = L_continuous, betaLong_linear = L_linear,
                     betaLong_aux = L_aux,
                     betaEvent_intercept = S_intercept, betaEvent_binary = S_binary,
                     betaEvent_continuous = S_continuous, betaEvent_assoc = S_assoc, betaEvent_aux = S_aux,
                     b_sd = b_sd, b_rho = b_rho, prob_Z1 = prob_Z1, mean_Z2 = mean_Z2, sd_Z2 = sd_Z2, 
                     max_yobs = yobs, max_fuptime = fuptime, balanced = bal,
                     family = gaussian, return_eta = FALSE,
                     seed = Seed, interval = c(1e-08, 200))

## simjm uses the provided simulation parameter values to create:
## dataSim[[1]] contains: time-fixed data
##                        id, Z1 (binary cov), Z2 (continuous cov), eventtime, status
## dataSim[[2]] contains: time-updated data
##                        id, Z1, Z2, eventtime, status, tij (meas time), Yij_1 (longitudinal cov)

## Keep the last Ntrain individuals as the test dataset
testSim <- list(dataSim[[1]][dataSim[[1]]$id > Ntrain*R,], 
                    dataSim[[2]][dataSim[[2]]$id > Ntrain*R,])
dataSim[[1]] <- dataSim[[1]][dataSim[[1]]$id <= Ntrain*R,]
dataSim[[2]] <- dataSim[[2]][dataSim[[2]]$id <= Ntrain*R,]

## Number each dataset so we can refer to them easily later
dataSim <- addDatasetNum(dat=dataSim, Ntrain=Ntrain) ##This helper fxn is in simSupportFxns.R
## test data will be one as there is only 1 test dataset
testSim[[2]]$dataSet <- 1


## ---- Transform the longitudinal covariate ---- ##

## keep a copy of the original longitudinal covariate for reference
dataSim[[2]]$Yij_orig <- dataSim[[2]]$Yij_1
testSim[[2]]$Yij_orig <- testSim[[2]]$Yij_1

## apply the transformation of your choosing
## the one shown here is "Transformation B" from the thesis
## but any transformation can be used here
dataSim[[2]]$Yij_1 <- (dataSim[[2]]$Yij_1 * 50*dataSim[[2]]$Z1)/(2*(dataSim[[2]]$Z2)) + 
                 (1-dataSim[[2]]$Z1)*(dataSim[[2]]$Yij_1)*2/(dataSim[[2]]$Z2)
testSim[[2]]$Yij_1 <- (testSim[[2]]$Yij_1 * 50*testSim[[2]]$Z1)/(2*(testSim[[2]]$Z2)) + 
                 (1-testSim[[2]]$Z1)*(testSim[[2]]$Yij_1)*2/(testSim[[2]]$Z2)


## ---- Create data for landmarking-style analysis ---- ##

testDatLM <- getLMdataset(datLong=testSim[[2]])
dataSimLM <- getLMdataset(datLong=dataSim[[2]])



## ---- Run models ---- #

## ----
## Fit a JM to the simulated training dataset and 
## predict on the validation data
##   Remember: for the joint model, you don't use the landmarking-style datasets or the LMEMs

modelJoint <- list()
for(i in 1:R){
  modelJoint[[i]] <- runJointM( i, trainDataS=dataSim[[1]][dataSim[[1]]$dataSet==i,],
                           trainDataL=dataSim[[2]][dataSim[[2]]$dataSet==i,],
                           testDataL =testSim[[2]], testDataS =testSim[[1]],
                           w=predHor, tseq=tseq)
}


## ----
## Fit the Cox Landmarking to the simulated training dataset and 
## predict on the validation data

modelCoxLM <- list()
for( i in 1:R){
  modelCoxLM[[i]] <- runCoxLM(indx=i, trainData <- filter(dataSimLM, dataSet==i),
                              testData=testDatLM, w=predHor, tseq=tseq)
}


## ----
## Fit the Super Learner Landmarking to the simulated training dataset and 
## predict on the validation data
modelSL_LM <- list()
for(i in 1:R){
  modelSL_LM[[i]] <- runSL_LM(indx=i, trainData <- filter(dataSimLM, dataSet==i), 
                              testData=testDatLM, w=predHor, tseq=tseq,
                              SL.lib=SL.libQuick2)
}


## Predicted survival probabilities corresponding to each of the R
## training datasets are stored as a list for each method.
## These can be compared using Brier score, C-index, etc.




