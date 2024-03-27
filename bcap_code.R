#####
options(mc.cores = parallel::detectCores()) 
library(rstan) 
library(MASS)
library(loo) 


mod = mod_bcap = rstan::stan_model("bcap_model.stan")  # compile the stan file  



### the main function for Bayesian CAP estimation 
bcap_estimation = function(data, 
                           mod=mod_bcap, ## this is a "stanmodel"
                           d=2, 
                           iter_warmup  = 700, 
                           iter_sampling= 2000, 
                           chains=1,
                           vb_approx=FALSE, 
                           vb_algorithm ="meanfield", vb_tol_rel_obj=0.005, 
                           vb_eta=0.1, vb_adapt_engaged=FALSE,
                           vb_iter=5000,  
                           show_messages=TRUE, 
                           init = NULL,  ## alternatively you can initialize using "cap"
                           warn = -1, 
                           lp_re=TRUE,   
                           CrI_probs=c(0.025,0.975))
{
  options(warn=warn)
  
  data$d = d
  data$lp_re = as.numeric(lp_re) 
  
  if(is.null(init)) init = "random" 
  
  cap.obj=NULL 
  if(init =="cap") {
    
    if((data$p>d+1)) {
      if(min(data$t_vec)<= data$p) {
        tmp = data$Y.list 
        for(i in 1:data$n) {
          
          n.aug = data$p-data$t_vec[i]+1  
          
          ## just for initialization of the chain, we will make the sample covariance matrix non-singular (by adding small noise)
          tmp[[i]] = rbind(tmp[[i]], matrix(rnorm(n.aug*data$p,0,0.01), n.aug, data$p)) 
        }
      }
      
      cap.obj = cap::capReg(Y=tmp, X=data$X.mat, nD=data$d, method=c("CAP"),CAP.OC=TRUE,ninitial=10)  
      
      init0 =  list(U =as.vector(cap.obj$gamma) + rnorm(length(as.vector(cap.obj$gamma)),0,0.01), 
                    beta = matrix(t(cap.obj$beta)[,-1], nrow=d),  
                    beta0 = as.array(log(diag(t(cap.obj$gamma)%*% data$whitening_mat%*%data$whitening_mat%*% cap.obj$gamma))), 
                    z = matrix(0, d, data$n),
                    L_Omega = diag(d),Omega_sd=as.array(rep(0.1,d))) 
      
      init = replicate(chains, init0, FALSE)
      if(vb_approx) init = init0 
      
    } else {
      
      init = 'random'
    
    }
  }
  
   
  ### run HMC or 'variational bayes' method to approximate the posterior  
  if(vb_approx){   
    
    ## variational bayes (default is meanfield ADVI)
    fit =  rstan::vb(mod, 
                     data = data, 
                     output_samples = iter_sampling, 
                     algorithm = vb_algorithm,
                     tol_rel_obj=vb_tol_rel_obj, iter=vb_iter, 
                     eta=vb_eta, adapt_engaged=vb_adapt_engaged,
                     init= init)
    
  } else { 
    
    ## run HMC 
    fit =  rstan::sampling(mod, 
                           data = data, 
                           chains = chains,
                           iter = iter_sampling, 
                           warmup = iter_warmup,
                           show_messages = show_messages,
                           init= init)
  } 
  
  if (vb_approx) {
    
    mean_accept = div =tree_hit =NULL
    
  } else {
    
    sampler_params = get_sampler_params(fit, inc_warmup = FALSE)
    mean_accept = mean(sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))) 
    div = sum(sapply(sampler_params, function(x) sum(x[, "divergent__"]))) 
    tree_hit = sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"]== 10))) 
  }
  
  list_of_draws = rstan::extract(fit)     
  
   
  
  ### extracted draws 
  draws_Gamma= draws_Gamma_tmp = list_of_draws$Gamma   
  draws_beta= draws_beta_tmp  = list_of_draws$beta
  draws_beta0= draws_beta0_tmp = list_of_draws$beta0
  draws_beta0_org= draws_beta0_org_tmp = list_of_draws$beta0_org
  draws_z = draws_z_tmp = list_of_draws$z   
  draws_Omega_sd = draws_Omega_sd_tmp = list_of_draws$Omega_sd   
  draws_Omega_cor = draws_Omega_cor_tmp = list_of_draws$Omega_cor    
  draws_Omega = draws_Omega_tmp = list_of_draws$Omega  
  draws_DfD = list_of_draws$DfD    
  draws_beta0_org = draws_beta0_org_tmp = list_of_draws$beta0_org  ## the objects re-mapped to the original space 
  
  ### matrix of log-likelihood for N records (column) and nsample samples (row).   
  log_lik = extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)  
  log_lik1 = extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)   
  log_lik_null = extract_log_lik(fit, parameter_name = "log_lik_null", merge_chains = TRUE)  
  
  ### waic    
  waic = loo::waic(log_lik1)$estimates[3,1] - loo::waic(log_lik_null)$estimates[3,1]  
   
  tmp = list_of_draws$Gamma   
  nsample = dim(tmp)[1]    ## number of posterior samples   
  n = nrow(data$X.mat)     ## number of units (subjects)  
  
  for (k in 1:d) {
    
    idx = which.max(abs(tmp[1,,k]))
    
    ## align the posterior samples based on the first post-warmup samples: we keep the sign of the element with the maximum absolute value throughout the posterior samples.
    for (s in 1:nsample) {
      tmp[s,,k] = tmp[s,,k]* as.numeric(sign(tmp[s,idx,k]))
    } 
    
  }
   
  lp_tmp = array(NA, dim = c(nsample, n, d))
  for (i in 1:n) {
    for (s in 1:nsample) {
      lp_tmp[s,i,] = cbind(draws_beta0_org_tmp[s,], draws_beta_tmp[s,,]) %*% data$X.mat[i,] + draws_z_tmp[s,,i]
    }
  }
  
  ### sort the components according to the between-subject variance explained 
  lp_m = apply(lp_tmp, c(2,3), mean)  # compute the posterior mean of log-variance 
  lp_v = apply(lp_m, 2, var)   # compute the between-subject variance of the mean log-variance 
  ind = order(lp_v, decreasing = TRUE) 
  var_lp = lp_v[ind] 
  var_lp_cumulated = cumsum(var_lp) 
  pve = var_lp_cumulated/max(var_lp_cumulated)  # proportion of the between-subject variance explained
  
  ### place the components in the "correct" order 
  for(k in 1:d) {
    draws_Gamma[,,k] = draws_Gamma_tmp[,,ind[k]]
    draws_beta[,k,] = draws_beta_tmp[,ind[k],]
    draws_beta0_org[,k] = draws_beta0_org_tmp[,ind[k]]
    draws_beta0[,k] = draws_beta0_tmp[,ind[k]]
    draws_Omega_sd[,k] = draws_Omega_sd_tmp[,ind[k]] 
    draws_z[,k,] = draws_z_tmp[,ind[k],]
    for(k2 in k:d) {
      draws_Omega[,k,k2] = draws_Omega_tmp[,ind[k],ind[k2]]
      draws_Omega_cor[,k,k2] = draws_Omega_cor_tmp[,ind[k],ind[k2]]
    }
  }
  
  ### posterior mean
  Gamma.hat = apply(draws_Gamma, c(2,3), mean) 
  beta.hat = apply(draws_beta, c(2,3), mean)  
  beta0.hat = apply(draws_beta0, 2, mean)  
  beta0_org.hat = apply(draws_beta0_org, 2, mean)  
  Omega_sd.hat = apply(draws_Omega_sd, 2, mean) 
  Omega.hat = apply(draws_Omega, c(2,3), mean) 
  Omega_cor.hat = apply(draws_Omega_cor, c(2,3), mean) 
  DfD.hat = mean(draws_DfD) 
  
  ### credible intervals 
  CrI_Gamma = CrI_beta = CrI_beta0_org =  CrI_Omega = NULL   
  for(k in 1:d) {
    CrI_Gamma = rbind(CrI_Gamma, t(apply(draws_Gamma[,,k], 2, 
                                         function(x) quantile(x, probs =CrI_probs))) )
    
    CrI_beta = rbind(CrI_beta, t(apply(as.matrix(draws_beta[,k,]), 2, 
                                       function(x) quantile(x, probs =CrI_probs))))
    
    CrI_beta0_org = rbind(CrI_beta0_org, t(apply(as.matrix(draws_beta0_org[,k]), 2, 
                                                 function(x) quantile(x, probs =CrI_probs))))
    
    for(k2 in 1:d)  CrI_Omega = rbind(CrI_Omega, quantile(draws_Omega[,k2,k], probs =CrI_probs))
  }
  
  
  list(fit=fit, 
       draws_Gamma = draws_Gamma,
       draws_beta = draws_beta,
       draws_z=draws_z,
       draws_Omega_sd = draws_Omega_sd,
       draws_Omega = draws_Omega,  
       draws_beta0_org=draws_beta0_org,
       draws_beta0=draws_beta0,
       draws_DfD = draws_DfD,  
       X.mat=data$X.mat,
       var_lp=var_lp, pve=pve, 
       chains = chains, 
       init = init, 
       Gamma.hat=Gamma.hat, 
       beta.hat=beta.hat, 
       beta0.hat=beta0.hat, 
       beta0_org.hat= beta0_org.hat, 
       Omega_sd.hat=Omega_sd.hat,
       Omega.hat=Omega.hat,  
       DfD.hat=DfD.hat, 
       CrI_Gamma = CrI_Gamma, 
       CrI_beta = CrI_beta,
       CrI_beta0_org =CrI_beta0_org, 
       CrI_Omega=CrI_Omega, 
       d=d, p=data$p, q=data$q, N=data$N, n=data$n, nsample=nsample,
       mean_accept = mean_accept, div=div, tree_hit=tree_hit,
       cap.obj=cap.obj, 
       waic=waic,   
       log_lik_null=log_lik_null,
       log_lik1=log_lik1,
       log_lik=log_lik, 
       ind=ind)
}



### data generation function used in the simulation examples in the paper
data.gen =  function(n=100, ## number of subjects 
                     nt=10, ## number of time points 
                     d=2,   ## number of signal components
                     p=10,  ## response dimension 
                     q=4,   ## number of covariates (excluding the intercept) 
                     sigma_z=0.5,  ## standard deviation of random effects z in the linear predictor
                     rho=0.4,      ## correlation in the random effects z in the linear predictor
                     sigma_noise=0.5,  ## standard deviation of noise component eigenvalues 
                     misspecified=FALSE)  
{
  
  A0 = rstiefel::rmf.matrix(matrix(rnorm(p*p),p,p))  
  GAMMA = A0[,1:d] 
  L = A0[,-c(1:d)] 
  L_i = array(dim = c(n,p,p-d))    ## subject-specific matrices 
  for (i in 1:n) { 
    noise.scale = diag(exp(sigma_noise*MASS::mvrnorm(1,rep(0,p-d),diag(rep(1,p-d)))/2))
    L_i[i,,] = L %*% rstiefel::rmf.matrix(matrix(rnorm((p-d)*(p-d)),p-d,p-d)) %*% noise.scale
  }
  
  GAMMA_i = array(dim = c(n,p,d))  ## subject-specific basis matrices  
  for (i in 1:n) {
    if (misspecified) {
      theta_i = runif(1, -pi/10, pi/10)
      GAMMA_i[i,,] = GAMMA %*% matrix(c(cos(theta_i), sin(theta_i), -sin(theta_i), cos(theta_i)), ncol=d)
    } else {
      GAMMA_i[i,,] = GAMMA
    }
  }
  
  beta0 = rep(0.1, d)   ## "intercept"  
  B = matrix(0, d, q)   ## regression coefficient matrix 
  B[1,1:q] = c(0.4, -0.5, 0.5, -0.5, rep(0.3, q-4))[1:q]  
  B[2,1:q] = c(-0.3, 0.4, -0.4, 0.4, rep(0.3, q-4))[1:q]   
  
  X = cbind(rbinom(n,1,0.5),matrix(rnorm(n*(q-1)),n,q-1))  # covariates  
  cov_names = paste0("X",1:q)
  colnames(X) = cov_names
  X_mat = cbind(X, id=1:n) %>% as.matrix  # excluding the intercept term 
  
  cor_z = diag(d) + matrix(rep(rho,d^2),d,d)-diag(rep(rho,d))
  sd_z = diag(rep(sigma_z, d))
  Omega = sd_z %*% cor_z %*% sd_z  
  z = MASS::mvrnorm(n, rep(0,d), Omega)  # subject-specific log-variance component 
  lp = beta0 + X %*% t(B) + z   # linear predictor (lp) for log-variance associated with the GAMMA component
  lp_var = apply(lp, 2, var)    
  
  Sigma = array(data=NA, dim= c(n, p, p))  
  for (i in 1:n) {
    Sigma[i,,] = GAMMA_i[i,,]%*%diag(exp(lp[i,]))%*%t(GAMMA_i[i,,]) + L_i[i,,]%*%t(L_i[i,,])
  }
  
  Y_tmp = Y_list = vector("list", length = n)
  for (i in 1:n) {
    Y_list[[i]] = scale(MASS::mvrnorm(nt, rep(0,p), Sigma[i,,]), 
                         center = TRUE, scale = FALSE)  ## mean-remove the outcomes for each subject 
    Y_tmp[[i]] = cbind(Y_list[[i]], id =i)
  }
  Y_mat = do.call(rbind, Y_tmp) %>% as.matrix 
  
  res_names = paste0("R",1:p)
  colnames(Y_mat)[1:p] = res_names
  t_vec = sapply(Y_list, nrow)  ## vector of number of time points
  
  ## Merge the response signals and subject info into one data frame. 
  dat = merge(Y_mat, X_mat, by = "id")
  
  S = array(data=NA, dim= c(n, p, p))     ## a list of sample covariance matrices 
  S_bar = matrix(data=0, nrow=p, ncol=p)  ## population level covarince matrix
  for(i in 1:n) {
    S[i,,] = var(Y_list[[i]])
    S_bar = S_bar+S[i,,]/n
  }
  
  ## S_bar will be used as a reference point for tangent-space parametrization
  svd_S_bar  = svd(S_bar)
  S_bar_sqrt = svd_S_bar$u %*% diag(sqrt(svd_S_bar$d)) %*% t(svd_S_bar$v)
  whitening_mat = svd_S_bar$u %*% diag(1/sqrt(svd_S_bar$d)) %*% t(svd_S_bar$v)
  
  
  list(B=B, GAMMA=GAMMA, L_i=L_i,  
       beta0=beta0, 
       sigma_z=sigma_z, 
       sigma_noise=sigma_noise, 
       Omega=Omega, 
       z=z, lp = lp, 
       dat=dat, Y_list=Y_list, X=X,
       cov_names=cov_names, res_names=res_names, 
       whitening_mat=whitening_mat,  ## the "whitening matrix" used in the tangent space parametrization
       S=S,          ## a list of individual-specific sample covariance matrices (used in the DfD criterion)
       S_bar=S_bar,  ## sample marginal covariance (to be used in the tangent space parametrization)
       S_bar_sqrt=S_bar_sqrt, 
       Sigma=Sigma,  ## a list of individual-specific true covariance matrices 
       p=p, d=d, q=q, N=nrow(dat), n=n, nt=nt,
       t_vec =t_vec, ## vector encoding the number of time points per subject
       lp_var=lp_var)
}



one_Experiment =  function(iter, 
                           d=2,   ## number of components 
                           n=100, ## number of subjects 
                           nt=10, ## number of time points
                           p=10,  ## response dimension 
                           q=4,   ## number of covariates (excluding the intercept) 
                           sigma_z=0.5,  ## standard deviation of random effects z in the linear predictor
                           rho=0.4,      ## correlation in the random effects z in the linear predictor 
                           sigma_noise=0.5, ## standard deviation of noise component eigenvalues 
                           misspecified= FALSE, 
                           vb_approx=FALSE, 
                           chains=1, 
                           iter_warmup= 700, 
                           iter_sampling= 2000,
                           lp_re=TRUE,
                           mod = mod_bcap, ## mod_bcap is a "stan" object   
                           beta0_sd=2.5,
                           beta_sd=2.5,
                           eta=1)
{
  
  set.seed(iter) 
  # generate data 
  dat.obj =  data.gen(n= n,
                      nt= nt,
                      p= p,
                      q= q, 
                      sigma_z=sigma_z, 
                      rho=rho,
                      sigma_noise=sigma_noise,
                      misspecified = misspecified)
  
  # organize the data for Stan 
  data =  list(p  = dat.obj$p,
               q  = dat.obj$q,
               N  = dat.obj$N, 
               n  = dat.obj$n, 
               id = dat.obj$dat[,"id"], # subject id
               X  = as.matrix(dat.obj$dat[,dat.obj$cov_names], ncol=dat.obj$q), 
               Y  = as.matrix(dat.obj$dat[,dat.obj$res_names], ncol=dat.obj$p),   
               whitening_mat = dat.obj$whitening_mat,    
               t_vec = dat.obj$t_vec, 
               S  = dat.obj$S,  
               Y.list = dat.obj$Y_list, 
               X.mat = cbind(1,dat.obj$X),  # for frequentist cap regression  
               beta0_sd =beta0_sd,
               beta_sd = beta_sd,
               eta = eta)
  
  # fit the model  
  time_elapsed = system.time({res = bcap_estimation(data=data, mod=mod, d=d,
                                                    chains=chains, 
                                                    iter_warmup = iter_warmup, 
                                                    iter_sampling=iter_sampling,
                                                    lp_re=lp_re)})[3]
  
  # some diagnostic measures  
  res$mean_accept
  res$div
  res$tree_hit
  
  list_of_draws   =  rstan::extract(res$fit) 
  draws_Gamma_tmp =  list_of_draws$Gamma
  draws_beta_tmp  =  list_of_draws$beta
  draws_mu_tmp  = list_of_draws$mu
  draws_sigma = list_of_draws$sigma
  draws_DfD = list_of_draws$DfD  
  
  if(d==dat.obj$d) {
    ## gamma rmse
    rmse.gamma1 =  1- abs(t(res$Gamma.hat[,1]) %*% dat.obj$GAMMA[,1]) 
    rmse.gamma2 =  1- abs(t(res$Gamma.hat[,2]) %*% dat.obj$GAMMA[,2]) 
    
    ## beta rmse 
    rmse.beta1 = sqrt(mean((res$beta.hat[1,] - dat.obj$B[1,])^2))
    rmse.beta2 = sqrt(mean((res$beta.hat[2,] - dat.obj$B[2,])^2))  
    
    # intercept rmse
    rmse.beta0 = sqrt(mean((res$beta0_org.hat - dat.obj$beta0)^2))
    
    ## Omega rmse   
    rmse.omega =sqrt(mean((as.vector(res$Omega.hat - dat.obj$Omega)[c(1,2,4)])^2))  
    
  } else {
    rmse.gamma1 = rmse.gamma2 = rmse.beta1 = rmse.beta2 = rmse.beta0 = rmse.omega = 0
  }
  
  
  ## coverage computation   
  coverage_GAMMA = rep(0, length(dat.obj$GAMMA))
  coverage_BETA  = rep(0, length(c(dat.obj$B[1,], dat.obj$B[2,])))
  coverage_beta0 = rep(0, length(dat.obj$beta0))
  coverage_Omega = rep(0, length(as.vector(dat.obj$Omega)[c(1,2,4)]))
  
  theta = dat.obj$GAMMA
  CrI = res$CrI_Gamma
  for(j in 1:nrow(CrI)){
    coverage_GAMMA[j] = 
      max(ifelse(CrI[j,1]<=  theta[j] &  theta[j]<= CrI[j,2], 1,0),
          ifelse(CrI[j,1]<= -theta[j] & -theta[j]<= CrI[j,2], 1,0))
  }
  
  theta = c(dat.obj$B[1,], dat.obj$B[2,])
  CrI = res$CrI_beta
  for(j in 1:nrow(CrI)){
    coverage_BETA[j] = ifelse(CrI[j,1]<=  theta[j] &  theta[j]<=CrI[j,2], 1,0)
  }
  
  coverage_gamma1 = mean(coverage_GAMMA[1:p])
  coverage_gamma2 = mean(coverage_GAMMA[1:p+(d-1)*p])
  coverage_beta1 = mean(coverage_BETA[1:q])
  coverage_beta2 = mean(coverage_BETA[1:q+(d-1)*q])
  
  if(d==dat.obj$d) {
    theta = dat.obj$beta0
    CrI = res$CrI_beta0_org
    for(j in 1:nrow(CrI)) {
      coverage_beta0[j] = ifelse(CrI[j,1]<=  theta[j] &  theta[j]<=CrI[j,2], 1,0)
    }
    
    theta = as.vector(dat.obj$Omega)[c(1,2,4)]
    CrI = res$CrI_Omega[c(1,2,4),]
    for(j in 1:nrow(CrI)) {
      coverage_Omega[j] = ifelse(CrI[j,1]<=  theta[j] &  theta[j]<=CrI[j,2], 1,0)
    } 
    coverage_omega = mean(coverage_Omega)
  }  
  
  
  return(c(
    rmse.gamma1=rmse.gamma1, 
    rmse.gamma2=rmse.gamma2, 
    rmse.omega=rmse.omega,
    rmse.beta1=rmse.beta1, 
    rmse.beta2=rmse.beta2, 
    rmse.beta0=rmse.beta0,  
    DfD.hat = res$DfD.hat,
    waic = res$waic,  
    coverage_GAMMA=coverage_GAMMA,  # each element-wise coverage
    coverage_BETA =coverage_BETA,
    coverage_beta0=coverage_beta0,
    coverage_Omega=coverage_Omega,
    coverage_gamma1=coverage_gamma1, # coverage averaged across the elements 
    coverage_gamma2=coverage_gamma2,
    coverage_beta1=coverage_beta1,
    coverage_beta2=coverage_beta2,
    coverage_omega=coverage_omega,
    time_elapsed= time_elapsed, 
    mean_accept = res$mean_accept,
    div = res$div, 
    tree_hit = res$tree_hit))
}

