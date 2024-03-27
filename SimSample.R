## Note: the most updated code is on Github: http://github.com/syhyunpark/bcap


######################
#####  Example  ######
######################

rm(list=ls(all=TRUE))
source('bcap_code.R')
 
n=100   ## number of subjects 
nt=10   ## number of time points
p=20    ## dimensionality of signal
sigma_z=0.5  ## standard deviation of random effects z in the linear predictor
rho=0.4      ## correlation in the random effects z in the linear predictor

### (see more specification in "data.gen")
### generate data
dat.obj  = data.gen(n= n,
                    nt= nt,
                    p= p, 
                    sigma_z=sigma_z, 
                    rho=rho)


# Y.list = dat.obj$Y_list
# S = array(data=NA, dim= c(n, p, p))     ## a list of sample covariance matrices 
# S_bar = matrix(data=0, nrow=p, ncol=p)  ## population level covarince matrix
# for(i in 1:n) {
#   S[i,,] = var(Y.list[[i]])
#   S_bar = S_bar+S[i,,]/n
# } 
### S_bar will be used as a reference point for tangent-space parametrization
# svd_S_bar  = svd(S_bar)
# S_bar_sqrt = svd_S_bar$u %*% diag(sqrt(svd_S_bar$d)) %*% t(svd_S_bar$v)
# whitening_mat = svd_S_bar$u %*% diag(1/sqrt(svd_S_bar$d)) %*% t(svd_S_bar$v)

### organize the data for Stan 
data =  list(p  = dat.obj$p,  ## number of regions 
             q  = dat.obj$q,  ## number of covariates (excluding intercept)
             n  = dat.obj$n,  ## number of subjects 
             N  = dat.obj$N,  ## number of total data points (including within-subject data)
             id = dat.obj$dat[,"id"], ## subject id
             t_vec = dat.obj$t_vec,   ## number of within-subject time points (length n vector, one number per subject)
             X  = as.matrix(dat.obj$dat[,dat.obj$cov_names], ncol=dat.obj$q), ## N x q design (covariate) matrix
             Y  = as.matrix(dat.obj$dat[,dat.obj$res_names], ncol=dat.obj$p), ## N x p response matrix 
             whitening_mat = dat.obj$whitening_mat, ## S_bar^{-1/2} (a matrix used to whiten the response signal)
             S  = dat.obj$S,  
             Y.list = dat.obj$Y_list,     ## for frequentist cap regression  
             X.mat = cbind(1,dat.obj$X),  ## for frequentist cap regression  
             beta0_sd= 2.5,  ## prior sd for intercept 
             beta_sd = 2.5,  ## prior sd for regression coefficent
             eta = 1)  ## LKJ hyperprior for correlations in random effects


## estimate the bayesian cap model 
res = bcap_estimation(data=data, mod=mod,  ## mod is the compiled "stanmodel" from 'bcap_code.R'. 
                      d=2, ## number of projection components 
                      iter_warmup = 200, ## iter_warmup  = 700
                      iter_sampling=500, ## iter_sampling=1300
                      init="cap")  ## initialized the chain with Zhao et al 2021's method; default is init="random"

res$Gamma.hat  ## posterior mean of p x d projection matrix 
res$beta.hat   ## posterior mean ofd x q regression coefficient matrix 
#res$CrI_beta  ## credible interval of the regression coefficients 
res$Omega.hat  ## posterior mean of the random effect d x d covariace matrix 
res$beta0_org.hat  ## posterior mean of the linear predictor's interecept vector 
ls(res)
res$waic  ## the waic value 

res$cap.obj ## Zhao et al (2021) cap regression object 
res$cap.obj$gamma




#####################################
#####  Simulation illustration ######
#####################################

scenarios =  expand.grid(n= c(100, 200, 300, 400),
                         nt= c(10, 20, 30),
                         p= c(10, 20),
                         misspecified=c(0,1))  # misspecified =1 means that there are no common covariate-relevant eigenvectors
scenarios
results.aggregated  = vector("list", length= nrow(scenarios))

set.seed(2023)
#for(scenario in 1:nrow(scenarios))
#{
scenario = 1
n = scenarios[scenario,1]
nt = scenarios[scenario,2]
p = scenarios[scenario,3] 
misspecified = scenarios[scenario,4] 

res = NULL 
n.rep=1 # 50

for(iter in 1:n.rep) {
  res_tmp = one_Experiment(iter=iter, 
                           n=n, nt=nt, p=p,
                           misspecified=misspecified)
  
  res = cbind(res, res_tmp)
  print(iter) 
  print("mean")
  print(round(apply(res, 1, mean),3))
  print("median")
  print(round(apply(res, 1, median),3))
  print("sd")
  print(round(apply(res, 1, sd),3))
}
#results.aggregated[[scenario]] = res

scenarios
print(scenario)
round(apply(res, 1, mean),3)
round(apply(res, 1, sd),3)
#}


######################################
########     End of code      ########
######################################