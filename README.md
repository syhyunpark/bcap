# R package for BCAP (Bayesian Estimation of Covariate Assisted Principal Regression)

The proposed Bayesian approach for estimating Covariate Assisted Principal (CAP) regression extends the frequentist method introduced by Zhao et al. (2021) by incorporating a probabilistic model coupled with a geometric formulation for dimension-reduced covariance objects. This framework allows for simultaneous inference on all relevant parameters, resulting in more interpretable results regarding the expression of covariate effects in the response vector's heteroscedasticity (covariance heterogeneity). Additionally, the dimension reduction approach applied to the outcomes eliminates the need to directly estimate subject-specific full p-by-p sample covariance matrices. This is advantageous in scenarios where the number of within-subject time points (volumes) is limited, as direct estimation can be prone to instability.


- Title: **Bayesian Estimation of Covariate Assisted Principal Regression**

- Authors: **Hyung Park<sup>a</sup> (parkh15@nyu.edu)**

- Affiliations:
   + 1. **Division of Biostatistics, Department of Population Health, NYU School of Medicine** 
  



## Setup Requirements
- R
- Stan (with rstan interface)


## Code Instructions

- The code for the proposed methodology is included in **bcap** folder. Please download the files in the folder to implement the method.
  + The main implementation functions are in **bcap_code.R**.
     + **bcap_model.stan** is a Stan file that specifies the model.
     + **SimSample.R** includes simulation examples demonstrating the implementation of the method.
     + **Application.R** contains code for applying the method and conducting data analysis procedures on a HCP (rs-fMRI PTN 820 subjects) dataset. 

### Main function: bcap_estimation
#### Arguments
+ `data`:  a **list** object, containing the data needed for Stan. This data is used as input for the Bayesian estimation process
   + `p`: dimensionality of the signal  
   + `q`: number of covariates (excluding intercept)
   + `n`: number of subjects 
   + `N`: number of total data points (including within-subject data)
   + `id`: vector of subject id (length N vector)
   + `t_vec`: number of within-subject time points (length n vector, one number per subject)
   + `X`: N x q design (covariate) matrix (excluding intercept)
   + `Y`: N x p response matrix 
   + `whitening_mat`: p x p matrix used to whiten the response signal (corresponding to $\bar{\Sigma}^{-\frac{1}{2}}$)
   + `S`: n x p x p list of subject-specific p x p covariance estimates; this is only to compute the Deviation from Diagonality (DfD) criterion and not essential for running BCAP 
   + `Y.list`: length-n list of $T_i$ x p observed signals; this is for frequentist cap regression (for initialization of the chain); otherwise set it to NULL
   + `X.mat`: n x (q+1) matrix of covariates (including the intercept); this is for frequentist cap regression (for initialization of the chain); otherwise set it to NULL
   + `beta0_sd`: prior sd for intercept
   + `beta_sd`: prior sd for regression coefficent
   + `eta`:  LKJ hyperprior for correlations in the random effects in the linear predictor
+ `mod`:  a **stanmodel** object compiled from the **bcap_model.stan** file
+ `d`: number of projection components (the rank of the dimension reduction matrix, Gamma)
+ `iter_warmup`: default is 700
+ `iter_sampling`: default is 2000 (which gives 1300 post-warmup samples)
+ `chains`: default is 1; can set up multiple Markov chains by specifying a different value for this parameter
+ `vb_approx`: default is FALSE; when set to TRUE, it conducts variational inference (additional adjustable parameters within the function **bcap_estimation** can further customize this process)
+ `init`: default is "random", meaning an MCMC chain initializes with a random value; alternatively, can be initialized using a fit from CAP regression by Zhao et al. by setting init to "cap"
+ `lp_re`: TRUE/FALSE indicator whether the subject-specific random effects are to be included in the linear predictor
+ `CrI_probs`: default is c(0.025, 0.975), representing the lower and upper bounds of a 95% credible interval
  
#### Value
+ a **list** object, containing draws from the MCMC-approximated posterior distributions of ``Gamma`` (dimension reduction matrix), ``beta`` (regression coefficient matrix), ``beta0`` (intercept, both in original parametrization and tangent parametrization), ``z`` (subject-level random effect), ``Omega`` (random effect covariance), ``DfD``, and the model waic-based expected deviance estimate value (``waic``), as well as several MCMC diagonostic statistics (such as ``mean_accept``, ``div``, ``tree_hit``) 

### The arguments of other functions are described within R files.

 
