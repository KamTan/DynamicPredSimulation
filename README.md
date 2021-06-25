# DynamicPredSimulation
Sample code illustrating a simulation study that compares three methods for dynamic survival prediction: joint modelling, Cox landmarking and Super Learner landmarking.

Code to accompany chapter 4 of the thesis by Kamaryn Tanner:
"Dynamic Prediction, Mediation and Communication for Survival Outcomes, with Applications to Cystic Fibrosis"
  
For questions and comments regarding the code, please contact Kamaryn Tanner (kamaryn.tanner1@lshtm.ac.uk). As code is updated, it will be posted to https://github.com/KamTan/. 
See also my repository github.com/KamTan/DynamicPrediction.

This code was written in R v4.0.2 running on a Windows 10 PC.

The code is organised in three files:

 * In simMain.R, we load the necessary libraries and source files and then simulate R training datasets, each with Ntrain individuals plus one test dataset using simjm from the simjm library.  Once the data has been created, each of the three methods is trained on each training dataset and then dynamic survival predictions are made on the test dataset.

 * simParams.R contains all of the parameters for controlling the simulation as well as settings for landmarking and the Super Learner.  By creating different versions of simParams.R, you can investigate different simulation scenarios without changing any code. 
 
 * simSupportFxns.R contains all of the helper functions called in simMain.R to implement the three methods.

For more information about dynamic survival prediction using a machine learning ensemble, please see 
"Dynamic Survival Prediction Combining Landmarking with a Machine Learning Ensemble: Methodology and Empirical Comparison" by Kamaryn T. Tanner, Linda D. Sharples, Rhian M. Daniel and Ruth H. Keogh. (2021) Journal of the Royal Statistical Society Series A, 184 (1) p3-30. https://doi.org/10.1111/rssa.12611


For more information about the simjm function used to generate simulated data from a joint model, please see:
https://github.com/sambrilleman/simjm


