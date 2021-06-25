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
##  This file contains the helper functions called from simMain.R
## ------------------------------------------------------------------
#####################################################################


## add a dataset number for use later
addDatasetNum<- function(dat, Ntrain){
  dat[[1]]$dataSet <- 0
  dat[[2]]$dataSet <- 0
  begRow <- 1
  begId <- 1
  for(i in 1:R){
    endId <- begId + Ntrain -1
    dat[[1]]$dataSet[begId:endId] <-i
    endRow <- max(which(dat[[2]]$id==endId)) 
    dat[[2]]$dataSet[begRow:endRow] <- i
    begId <- endId +1
    begRow <- endRow +1
    
  }
  return(dat)
}

## get the datasets into a list to make them easier to work with
dataSimToList <- function( dat){
  dL <- list()
  for( i in 1:R){
    dL[[i]] <- dat[dat$dataSet==i,]
  }
  return( dL)
}


## This function takes a dataset with Yij at tij and xforms it to one with 
##  tij at LM=0, 1, 2, etc. based on LMgridSize
getLMdataset <- function( datLong){
  datBase <- ddply(datLong,.(id), tail,1)
  l <- which(datBase$eventtime < LMfirst)
  if(length(l)>0){
    datBase <- datBase[-l,]
  }
  ##pre-allocate
  ids <- unique(datBase$id)
  dim <- length(ids)*(LMlast-LMfirst+1)*(1/LMgridSize)
  datLM <- data.frame(id=rep(-99,dim), eventtime=rep(-99,dim), status=rep(-99,dim),
                      Z1=rep(-99,dim), Z2=rep(-99,dim), tij=rep(-99,dim), dataSet=rep(-99,dim))
  ##track your row
  rowBeg <- 1
  rowEnd <- 1
  for(i in ids){
    maxLM <- min( LMlast, floor(datBase$eventtime[datBase$id==i]))
    rowEnd <- rowBeg + (maxLM -LMfirst)    
    datLM$id[rowBeg:rowEnd] <- i
    datLM$eventtime[rowBeg:rowEnd] <- datBase$eventtime[datBase$id==i]
    datLM$status[rowBeg:rowEnd] <- datBase$status[datBase$id==i]
    datLM$Z1[rowBeg:rowEnd] <- datBase$Z1[datBase$id==i]
    datLM$Z2[rowBeg:rowEnd] <- datBase$Z2[datBase$id==i]
    datLM$tij[rowBeg:rowEnd] <- seq(LMfirst,maxLM)
    datLM$dataSet[rowBeg:rowEnd] <- datBase$dataSet[datBase$id==i]
    rowBeg <- rowEnd + 1
  }
  ##now get rid of any extra rows
  datLM2 <- datLM[1:rowEnd,]
  
  ##We also have to add in the "observed" Yij and its measurement time 
  myCols=c("id", "eventtime", "status", "Z1", "Z2", "tij", "Yij_1", "Yij_orig")
  datLM3 <- full_join(datLM2, datLong[,myCols], by=c("id", "eventtime", "status", "Z1", "Z2", "tij"))
  
  ##Sort it by id then by timestamp
  datLM3 <- datLM3 %>% arrange(id, tij)
  
  Yij_1_locf <- by(datLM3, datLM3$id, function(df){na.locf(df$Yij_1, na.rm=F)})
  datLM3$Yij_1_locf <- unlist(Yij_1_locf)
  
  return(datLM3)
}

#######################
#  For joint model
#######################

## Note this function doesn't take w because I've had to hard-code it at 2 and 5
runJointM<-function(indx, trainDataS, trainDataL, testDataL, testDataS, w, tseq){
  
  ids <- unique(trainDataL$id)
  trainDataS <- trainDataS[trainDataS$id %in% ids,]
  #Fit model for the longitudinal process
  lmeFit <- lme(Yij_1 ~ tij, random =~tij | id, data=trainDataL,
                control=lmeControl(opt='optim', msMaxIter=100))
  #Then the survival model
  coxFit <- coxph(Surv(time=eventtime, event=status) ~ Z1 + Z2,
                  data = trainDataS, x=TRUE, model=TRUE)
  # Third, the joint model
  #  Using Weibull because it is very fast but may need to consider other
  #  specifications such as method="spline-PH-GH"
  jointFit <- jointModel(lmeFit, coxFit, timeVar="tij")
  
  # Fourth, get predictions at same grid of landmark times used in other methods
  cat("\nGetting joint model predictions  ") #Provide info on progress
  cat(indx)
  endAge <- endAgeJM
  testIds <- unique(testDataL$id)
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds

  for( j in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    ## limit to ids that were at risk and only predict using data before landmark time
    sData <- testDataL[testDataL$eventtime>tseq[j],]
    sData <- sData[sData$tij <= tseq[j],]
    sIds <- as.character(unique(sData$id))
    predProbs <- matrix(nrow=3, ncol=length(sIds))
    colnames(predProbs) <- sIds
    #This loop calculates the conditional survival probabilities by id
    for( l in 1:length(sIds)){
      ND <- sData[sData$id %in% sIds[l],]
      ND$LM <- tseq[j]
      p1 <- JM::survfitJM(jointFit, newdata=ND, survTimes=c(tseq[j], (tseq[j]+w)), last.time = "LM")
      ##capture mean survival probability (column 2)
      sp <- p1[[1]][[1]][,2]
      if(p1[[1]][[1]][1,1] != tseq[[j]]){ #then we need to add a 1 to the front of sp
        sp <- c(1.0, sp)
      }
      testSurvPreds[j, which(testIds==sIds[l])] <- t(sp)[2]
    }
  }
  ##Predicted survival probabilities on the test set for this training dataset
  return(testSurvPreds)  
}

  
runCoxLM<-function(indx, trainData, testData, w, tseq){
  
  ## get stacked data for landmarking
  lmData <- genStackedData(dat = subset(trainData, select=-c(dataSet, Yij_1,Yij_orig)), 
                           w, tseq)
  ##run Cox regression stratified on landmark time
  LMcox <- coxph(Surv(lmTime, eventtime, status) ~ Z1 + Z2 + Yij_1_locf + 
                       strata(lmTime) + cluster(id),data=lmData, method="breslow")

  ##Get predictions along landmark time grid
  testIds <- unique(testData$id)
  ##set up storage for predicted survival probabilities
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds
  
  for( k in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    sData <- testData[testData$eventtime>tseq[k],]
    sData <- sData[sData$tij <= tseq[k],]
    sData$eventtime[sData$eventtime > (tseq[k]+w)] <- (tseq[k]+w)
    sData <- ddply(sData,.(id), tail,1)  # Keep only the most recent record for each id
    newData <- subset(sData, select= c(Z1, Z2, Yij_1_locf))
    
    ## Predicted curves from a coxph model have one row for each stratum in the Cox model fit 
    sfhor <- survfit(LMcox,newdata=newData)
    ##The data I want will be between indices:
    if(k == 1){
      index1 <- 1  
    }else {
      index1 <- cumsum(sfhor$strata)[k-1] + 1
    }
    index2 <- index1 + sfhor$strata[k] - 1
    
    nhor <- index2
    ids <- unique(sData$id)
    if( length(ids)>1){
      for(j in 1:length(ids)){#loop over individuals
        testSurvPreds[k, which(testIds==ids[j])] <- sfhor$surv[nhor,j]
      }
    }else {
      testSurvPreds[k,ids[1]] <- sfhor$surv[nhor]
    }
  }
  ##Predicted survival probabilities on the test set for this training dataset
  return(testSurvPreds)    
}

runSL_LM <- function(indx, trainData, testData, w, tseq, SL.lib){
  
  ## get stacked data for landmarking with SL
  sd <- genStackedSLData(dat=subset(trainData, select=-c(dataSet,Yij_1,Yij_orig)), w, tseq)
  stackedData.disc <- sd
  
  SLmodel <- stackedSL(SL.lib, stackedData= stackedData.disc,  
                       verbose=verbose, numFolds=numSLFolds, 
                       method=SL.method, indx=indx)
  
  ##Get predictions along landmark time grid
  testIds <- unique(testData$id)
  ##set up storage for predicted survival probabilities
  testSurvPreds = matrix(NA, length(tseq), length(testIds))
  colnames(testSurvPreds) <- testIds
  colN <- paste0("LM", tseq)
  colN2 <- paste0("LM", tseq, "sq")
  
  for( k in 1:length(tseq)){ #loop over grid of landmark times to make predictions at
    
    thor <- tseq[k]+w
    thisLM <- paste0("LM", tseq[k])
    thisLM2 <- paste0("LM", tseq[k], "sq")
    sData <- testData[testData$eventtime>tseq[k],]
    sData <- sData[sData$tij <= tseq[k],] # delete all records where measurement collected after tseq[i]
    sData <- ddply(sData,.(id), tail,1) # Keep only the most recent record for each id
    ids2 <- as.character(unique(sData$id))
    delta.bound <- sort(unique(stackedData.disc[thisLM][stackedData.disc[thisLM]!=0] ))
    nhor <- length(delta.bound)
    ## Must put the data for prediction into the same format as we used to fit the model
    ND <- createPredictData(dataX= subset(sData[sData$id %in% ids2,], 
                                          select= c(Z1,Z2,Yij_1_locf)),time=delta.bound)
    
    ## We need to know delta.bound but don't want it in the data
    ND2 <- ND[1:(ncol(ND)-1)]
    ## add the landmark time dummy variables
    dum <- matrix(0, nrow=nrow(ND2), ncol=length(tseq))
    colnames(dum) <- colN
    rownames(ND2) <- NULL
    ND2 <- cbind(ND2, dum)
    # Add the landmark time as an interaction with time
    ND2[,thisLM ]<- ND$time   
    sfhor <- predict(SLmodel, newdata = ND2 ) 
    ## sfhor contains pred and library.predict -- these are conditional hazards
    ## there will be one row per row in dhor.disc
    #loop over individuals to calc survival probs
    for(j in 1:length(ids2)){
      beg <- ((j-1)*nhor) +1
      en <- j*nhor
      survP <- cumprod((1-sfhor$pred[beg:en]))
      testSurvPreds[k, ids2[j]] <- survP[nhor]
    }  
  }
  return(testSurvPreds)
}

  
  
 
genStackedData <- function(dat, w, tseq){
  
  stackedData <- data.frame()
  for(i in 1:length(tseq)){
    dat_i <- dat[dat$eventtime>tseq[i],]   #keep only those still at risk at s  
    dat_i$status[dat_i$eventtime>(tseq[i]+w)] <- 0   #ignore events after s+w
    dat_i$eventtime[dat_i$eventtime>(tseq[i]+w)] <- (tseq[i]+w)  #administratively censor at s+w
    dat_i <- dat_i[dat_i$tij <= tseq[i],] # delete all records where measurement collected after s
    dat_i <- ddply(dat_i,.(id), tail,1)  # Keep only the most recent record for each id
    dat_i$lmTime <- tseq[i]                   # Add a column for the landmark time
    stackedData <- rbind(stackedData, dat_i) 
  }
  return(stackedData)
}   ##closes function

##
## --- Function to create stacked dataset for Super Learner landmarking
##     The difference is that this dataset will be discretised
##

genStackedSLData <- function(dat, w, tseq){
  n.delta <- 5  # number discrete intervals to use
  stackedData <- data.frame()
  colN <- paste0("LM", tseq)
  for (i in 1:length(tseq)) {
    dat_i <- dat[dat$eventtime>tseq[i],]   #keep only those still at risk at s
    dat_i$status[dat_i$eventtime>(tseq[i]+w)] <- 0   #ignore events after s+w
    dat_i$eventtime[dat_i$eventtime>(tseq[i]+w)] <- (tseq[i]+w)  #administratively censor at s+w
    dat_i <- dat_i[dat_i$tij <= tseq[i],] # delete all records where measurement collected after s
    dat_i <- ddply(dat_i,.(id), tail,1) #get the last record only
    
    # Discretize it based on quantiles of the event times
    delta.upper <- createDiscreteIntervals(time= dat_i$eventtime, event=dat_i$status, 
                                           n.delta=n.delta)
    dat_i.newX <- createDiscrete2(time= dat_i$eventtime, event=dat_i$status, 
                                  dataX = subset(dat_i, select=-c(eventtime, status, tij)), 
                                  delta.upper = delta.upper, s=tseq[i])
    
    dat_i.disc <- dat_i.newX[!is.na(dat_i.newX$N.delta), ] 
    
    ## Add dummies for landmark time and interact with delta.lower
    dum <- matrix(0, nrow=nrow(dat_i.disc), ncol=length(tseq))
    colnames(dum) <- colN
    
    dat_i.disc <- cbind(dat_i.disc, as.data.frame(dum))
    thisLM <- paste0("LM", tseq[i])
    dat_i.disc[,thisLM ]<- dat_i.disc$delta.lower   # Add the landmark time as an interaction with time
    
    #rbind is known to be slow...
    stackedData <- rbind(stackedData, dat_i.disc)
  } 
  
  return(stackedData)
  
}   

##
## --- Function to determine the boundaries of each discrete interval
##

createDiscreteIntervals <- function(time, event, n.delta) {
  # time is the vector of survival times
  # n.delta is the number of discrete intervals
  n <- length(time)
  n.delta <- min(n.delta, length(unique(time[event==1])))
  probs.delta <- seq(from=0, to=1, length.out=(n.delta+1))
  delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
  return(delta.upper[-1])	
}

## This function was adapted from code from the Github account of Eric Polley
## https://github.com/ecpolley/SuperLearner_Old/blob/master/R/createDiscrete.R
## Also see Chapter 16 "Super Learning for Right-Censored Data" in the book
## "Targeted Learning: Causal Inference for Observational and Experimental Data"
## by van der Laan and Rose, Springer, 2011.

createDiscrete2 <- function(time, event, dataX, delta.upper, s) {
  
  n <- length(time)
  delta.lower <- c(s, delta.upper[-length(delta.upper)])
  n.delta <- length(delta.upper)
  ID <- dataX$id
  dat <- cbind(ID, time, event, dataX)
  
  ## Note: Can use an apply statement here instead of pre-allocating and looping but
  ## apply converts data.frame to a matrix and uses least restrictive class which may be character
  long.dat <- matrix(NA, n.delta*n, ncol(dat))
  long.dat <- as.data.frame(long.dat)
  colnames(long.dat) <- colnames(dat)
  for( k in 1:nrow(dat)){ ##loop over ids
    b <- (n.delta*(k-1))+1
    e <- n.delta*k
    long.dat[b:e,] <- dat[k,]
  }
  
  N.delta <- rep(NA, nrow(long.dat))
  long.dat <- cbind(long.dat, delta.lower, delta.upper, N.delta)
  # Include censored people in the interval in which they were censored  
  long.dat$N.delta <- ifelse(long.dat$time > long.dat$delta.upper, 0, 
                             ifelse(long.dat$event==1, ifelse(long.dat$time <= long.dat$delta.lower, NA, 1), 
                                    ifelse(long.dat$time>long.dat$delta.lower, 0, NA)))
  
  m <- delta.upper[n.delta]
  long.dat$N.delta <- ifelse(long.dat$time == m & long.dat$delta.upper == m, 
                             ifelse(is.na(long.dat$N.delta), 0, long.dat$N.delta), long.dat$N.delta)
  
  return(long.dat)	
}


##
## --- This function is similar to the above but formats data we will make predictions on
##

createPredictData <- function(dataX, time){
  # Need to create a prediction set with the relevant data for this window
  n <- nrow(dataX)
  n.delta <- length(time)
  
  long.DATA <- apply(dataX, 2, function(x) rep(x, times=rep(n.delta, times=n)))
  long.DATA <- cbind(long.DATA, time)
  long.DATA <- as.data.frame(long.DATA)
  return(long.DATA)	
} 

##
## --- Function to apply SuperLearner to our prepared dicretised stacked dataset
##

stackedSL <- function(SL.lib, stackedData, newXdata=NULL, verbose=FALSE, 
                      numFolds=10, method="method.NNLS", indx){  
  
  
  ## fit the SuperLearner to the single stacked dataset
  set.seed(12)
  df.X <- subset(stackedData, select= -c(ID, time, event, id, delta.lower, delta.upper, N.delta ))
  
  cat("\nBegin SuperLearner ")
  cat(indx)
  cvControl = list(V=numFolds)
  
  SL.hor<- SuperLearner( Y = stackedData$N.delta, X = df.X,  
                         SL.library=SL.lib, 
                         method=method, verbose=verbose,
                         id=stackedData$id, family=binomial(), cvControl=cvControl)
  
  return(SL.hor)
}
