# R package for BCAP (Bayesian Estimation of Covariate Assisted Principal Regression)

Extending the frequentist approach developed in Zhao and others (2021) under a probabilistic model (coupled with a geometric formulation of the dimension-reduced covariance objects), the proposed Bayesian estimation approach to Covariate Assisted Principal (CAP) regression provides a framework to conduct inference on all relevant parameters simultaneously, that produces more interpretable results regarding how the covariates’ effects are expressed in the response vector. Furthermore, the outcome dimension reduction approach avoids the need to work with subject-specific full p-by-p sample covariance matrices, which can suffer from estimation instability when the number of time points (volumes) is not large (which is typically the case for fMRI signals). Generally, the CAP formulation allows for a more targeted and efficient analysis by identifying the specific components of the outcome data relevant to the association between covariates and functional connectivity.


- Title: **Bayesian Estimation of Covariate Assisted Principal Regression**

- Authors: **Hyung Park<sup>a</sup> (parkh15@nyu.edu)**

- Affiliations:
   + 1. **Division of Biostatistics, Department of Population Health, NYU School of Medicine** 
  



## Setup Requirements
- R
- Stan

<!--- 
## Code Instructions

- The code for the proposed methodology is included in **ICATemporalNetwork** folder. Please download all the files in the folder to implement the method.
  + The main function for the method is **ICATemporalNet.R** and **ICATemporalNetBoot.R** which allows bootstraps.
  + To use **ICATemporalNet.R**, it requires **estICA.R** which estimates the non-Gaussian signals and then removes them from raw data, and **temporalNet.R** which estimates the temporal network. 


 
- **Examples** folder contains examples.
   + **genData.R**: generate simulated data
   + **example.R**: an example to implement the method
   + **Sim_Scenario1.R**: simulations of Scenario 1
   + **Sim_Scenario2.R**: simulations of Scenario 2

### Main function: ICATemporalNet
#### Arguments
+ `Yts`: input data, the user must supply a list of Yts, where each element is a N*K data matrix at time t. N is sample size, K is the number of nodes.
+ `N`: sample size
+ `Ntime`: total number of time points
+ `ncomp`:  maximum number of independent components to be chosen
+  `Ta`: use t<=Ta time points to estimate temporal network A
+  `Tc`: ues t>Tc time points to estimate contemporaneous network Gamma

#### Value
+ `estIC`: results from non-Gaussian estimation step. output from estICA.R
+ `estRts`: R(t),residuals after removing non-Gaussian signals
+ `estS`: independent components S
+ `estUts`: non-Gaussian signals U(t)
+ `estWt`: weight matrix w(t)
+  `nIC`: number of independent components
+  `A`: temporal network
+  `Gamma`: contemporaneous network
+  `Omega`: covariance matrix of e(t), inverse of Gamma

### The arguments of other functions are described within R files.
-->
 
