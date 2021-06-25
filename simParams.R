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
##  This file contains parameters and settings used to control the generation
##   of simulated data as well as certain settings for landmarking and
##   Super Learner.  By creating different versions of this file, different
##   simulation scenarios can be tested.
## ------------------------------------------------------------------
#####################################################################

## ----
## Set parameters to control the simulation

# used for naming files
fnSim <- "A.rds"

R <- 5               # How many simulated datasets?
Ntrain <- 100        # Number of subjects in training data & test data
Seed <- 1212

# Covariates
prob_Z1 <- 0.5    #binary covariate
mean_Z2 <- 50     #continuous covariate, normally distributed
sd_Z2 <- 1

b_sd = c(44, 3)  #The vector of standard deviations for the individual-level random effects.
b_rho = -0.9         #Correlation between the individual-level random effects. only relevant when >1 random effect

# Longitudinal submodel parameters
L_intercept <- 60
L_binary <- 0       #True coefficient for binary covariate
L_continuous <- 0    #True coefficient for continuous covariate
L_linear <- -3    #True coefficient for fixed effect linear term
L_aux <- 6.6         #True sigma for Gaussian model

# Survival submodel parameters
S_intercept <- -7.0  
S_binary <- 0.5     #True log hazard ratio for the binary covariate in the event submodel
S_continuous <- 0.07  #True log hazard ratio for the continuous covariate
S_assoc <- -0.15       #True log hazard ratio for the longitudinal process
S_aux <- 2.0     #True parameter value for shape parameter for Weibull model

# Details of longitudinal data to be simulated
yobs <- 20        #Maximum allowed number of longitudinal measurements for the biomarker 
fuptime <- 16     #The maximum follow up time. This time will also be used as the censoring time 
bal	<- FALSE      #FALSE means meas times will be chosen randomly from a uniform distribution [0,max_fuptime] 

# 
nameSurvTime = "eventtime"
nameEventInd="status"
nameMeasDate="tij"

# Settings for landmarking
predHor <- 5        # The prediction horizon(s)
endAgeJM <- fuptime    # The end age we use for getting joint model prediction
#grid of landmark times used for methods and evaluation
tseq <- seq(5,10,by=1)
LMfirst <- 5
LMlast <- 10
LMgridSize <- 1

# Things for Super Learner
# Note: this library is a minimal set of algorithms designed to run quickly
#       better performance may be achieved by adding other algorithm/tuning parameter
#       combinations
SL.libQuick2 <- list("SL.glm", "SL.gam", "SL.ranger", "SL.xgboost", "SL.mean")
numSLFolds <- 5  ##Number of x-validation folds in Super Learner. 10 is preferred
SL.method <- "method.NNLS"  ## implies a squared error loss function
verbose=FALSE
