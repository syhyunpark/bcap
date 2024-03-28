### Bayes CAP code for Application    
### the whole analysis and plots 
 
library(ggplot2)
library(reshape2)
library(cowplot)
library(magrittr)
library(stringr)
library(dplyr) 
library(rstan)
library(abind)
library(png)  
library(coda) 
library(data.table)
library(bayesplot)
library(gridExtra)
library(grid) 
source("bcap_code.R") ## allow some time, approximately 30 seconds, for compiling the Stan file

## Load HCP Data

###########################################################################################################################################
## The dataset (located in the folder named "HCP_PTN820") and some code snippets for reading in the HCP data are sourced from 
## https://github.com/ChristofSeiler/CovRegFC_HCP (Christof Seiler and Susan Holmes). 
## Please note that the dataset is provided in a compressed format and 
## needs to be extracted from the zip file named "HCP_PTN820.zip" before use.
###########################################################################################################################################

## Set appropriate working directory  
# setwd("/Users/hyung/Dropbox/BayesCAP/BCAP_Biostatistics_revision/code_Github")
subject_info = read.csv("HCP_PTN820/sample_info.csv")
colnames(subject_info)
dim(subject_info)
subject_info = subject_info[,c(
  "Subject",
  "Age",
  "Gender",
  "Acquisition",
  "PSQI_AmtSleep", # amount of sleep in hours
  "PSQI_Score" # Pittsburgh Sleep Quality Index (PSQI) Completed
)] 

# * short sleepers: average equal or less than 6 hours each night
# * conventional sleepers: average between 7 and 9 hours each night
sleep_duration = rep("undefined", nrow(subject_info))
sleep_duration[subject_info$PSQI_AmtSleep <= 6] = "short"
sleep_duration[(subject_info$PSQI_AmtSleep >= 7) & (subject_info$PSQI_AmtSleep <= 9)] = "conventional"
subject_info$sleep_duration = factor(sleep_duration)
table(subject_info$sleep_duration)
str(subject_info)
dim(subject_info)

table(subject_info$Age)
tmp = factor(subject_info$Age, levels =c("22-25","26-30","31-35","36+"),
              labels= c("22", "26", "31", "36"))
subject_info$age =  scale(as.numeric(as.character(tmp)))


## Load timeseries.
num_regions = 15
channel_names = paste0("R",1:num_regions)
path = paste0("HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d",num_regions,"_ts2")
file_names = list.files(path = path, pattern = ".txt")
file_subject_ids = strsplit(file_names, split = ".txt") %>% unlist
ts = lapply(file_subject_ids, function(Subject) {
  print(paste("reading subject:", Subject))
  full_path = paste0(path, "/", Subject, ".txt")
  timeseries = read.csv(full_path, header = FALSE, sep = " ")
  timeseries$Subject = Subject
  timeseries$run = lapply(1:4, function(i) rep(i,nrow(timeseries)/4)) %>% unlist
  timeseries
}) %>% do.call(rbind,.) %>% data.frame
names(ts)[1:num_regions] = channel_names 
dim(ts)

## Merge timeseries and subject info into one data frame.
ts_subject_info = merge(ts, subject_info, by = "Subject")
dim(ts_subject_info)


## Plot timeseries for one subject.
subject_names = names(table(ts_subject_info$Subject))
## for one typical subject, extract "timeseries"
timeseries = subset(ts_subject_info, Subject == subject_names[1])
## each subject has 1200 time points per session; indicate 1:1200 time points per session
timeseries$timepoint = rep(1:(nrow(timeseries)/4), 4)
timeseries$run = lapply(paste0("run_",1:4), function(i) rep(i,nrow(timeseries)/4)) %>% unlist

timeseries_long = reshape2::melt(timeseries, id.vars = c("timepoint","run", names(subject_info)))
ggplot(timeseries_long, aes(x = timepoint, y = value, color = variable)) + 
  geom_line() + facet_wrap(~run)
timeseries_long_subset = subset(timeseries_long,
                                timepoint < 51 & (variable == "R1" | variable == "R2" | variable == "R3"))
ggplot(timeseries_long_subset,aes(x = timepoint,y = value,color = variable)) + 
  geom_line() + facet_wrap(~run)



## What is the effective samples size: computed using autocorrelation function.
## (this takes quite some time, several minutes) 
## (if wish to skip this step, go right to: tp_per_subject = c(34,34,34,34) )
run = cut_number(1:nrow(timeseries), 4, labels=FALSE)  # length 4800
ess = lapply(subject_names, function(subject_name) {
  sapply(1:4, function(i) {
    timeseries = subset(ts_subject_info, Subject == subject_name)
    Y = timeseries[run==i,channel_names]
    coda::effectiveSize(Y) %>% min
  })  #%>% sum
}) #%>% unlist
ess_session = matrix(ess %>% unlist, ncol=4, byrow=TRUE) 
tp_per_subject = round(apply(ess_session, 2, min)) 
#tp_per_subject = c(34,21,25,31)  

tp_per_subject = c(34,34,34,34)
dim(ts_subject_info)
dim(timeseries)
colnames(ts_subject_info)

## Subsample time points within each subject to account for dependencies.
## We will subsample and store the data for each session separately.
ts_subject_info_subset_4_sessions = vector("list", length=4)
run = unlist(lapply(1:4, function(x) rep(x, tp_per_subject[x]))) 
for(run_id in 1:4)
{
  ts_subject_info_subset = subset(ts_subject_info,
                                  run==run_id & (sleep_duration == "short" | sleep_duration == "conventional")) %>% droplevels
  subject_names = names(table(ts_subject_info_subset$Subject))
  subsample_ids = lapply(subject_names, function(subject_name) {
    ids = which(ts_subject_info_subset$Subject == subject_name)
    seq(from = min(ids), to = max(ids), length.out = tp_per_subject[run_id]) %>% floor
  }) %>% unlist 
  ts_subject_info_subset_4_sessions[[run_id]] = ts_subject_info_subset[subsample_ids,] %>% droplevels
}

ts_subject_info_subset = ts_subject_info_subset_4_sessions[[1]]
dim(ts_subject_info_subset)

dim(ts_subject_info_subset_4_sessions[[1]])
#dim(ts_subject_info_subset_4_sessions[[2]])
#dim(ts_subject_info_subset_4_sessions[[3]])
#dim(ts_subject_info_subset_4_sessions[[4]])


######################
### HCP Session 1. ###
######################
data = ts_subject_info_subset_4_sessions[[1]]   
dim(data)
#[1] 24820    24
region_ids = str_detect(names(ts_subject_info),"R") %>% which
subject_ids = table(data$Subject) %>% names
n = length(subject_ids)
p = length(region_ids)

## collect subject-specific information 
dat = list()
dat$S = dat$corr = array(dim = c(n,p,p)) 
dat$logvar.diag = array(dim=c(n,p))
dat$corr.lower.tri = array(dim=c(n,p*(p-1)/2))
dat$x = dat$Yc = vector("list", length=n)

for(i in 1:n) {
  dat$x[[i]] = subset(data, Subject == subject_ids[i])[1,c("sleep_duration","Gender")]
  dat$Yc[[i]] = Yc = scale(subset(data, Subject == subject_ids[i])[,region_ids], 
                           center = TRUE, scale = FALSE)  # center the subject-specific brain signals
  dat$S[i,,] = (t(Yc)%*%Yc)/nrow(Yc) 
  dat$corr[i,,] = cor(Yc)   # only to conduct edge-wise regression  
  dat$corr.lower.tri[i,] = dat$corr[i,,][lower.tri(dat$corr[i,,])]
  dat$logvar.diag[i,] = log(diag(dat$S[i,,]))
}

dim(dat$corr.lower.tri)
#[1] 730 105 
FisherZ = function(rho){ 0.5*log((1+rho)/(1-rho)) } 
FZ.corr = apply(dat$corr.lower.tri, MARGIN =c(1,2), FisherZ)   
logvar.diag = dat$logvar.diag 

y.mat = cbind(FZ.corr, logvar.diag)
X.mat = model.matrix(as.formula("~sleep_duration + Gender + sleep_duration:Gender"), 
                     data=as.data.frame(do.call(rbind, dat$x))) 

## write a contrast vector 
K1 = matrix(c(0, 1, 0, 1), 1) #  1) short vs. conventional sleeper among male
K2 = matrix(c(0, 1, 0, 0), 1) #  2) short vs. conventional sleeper among female;
K3 = matrix(c(0, 0, 1, 1), 1) #  3) male vs. female among short sleeper; and
K4 = matrix(c(0, 0, 1, 0), 1) #  4) male vs. female among conventional sleeper.

##########################################
## first, run a element-wise regression ##
##########################################
library(multcomp)
p_values.j = estimates.j = rep(NA,4) 
p_values = estimates = NULL 
for(j in 1:ncol(y.mat)){
  dat.j = data.frame(y=y.mat[,j], X.mat[,-1]) 
  fit.j = lm(y~.,data=dat.j) 
  glht1 = summary(glht(fit.j, linfct = K1))$test 
  glht2 = summary(glht(fit.j, linfct = K2))$test 
  glht3 = summary(glht(fit.j, linfct = K3))$test 
  glht4 = summary(glht(fit.j, linfct = K4))$test 
  estimates.j[1] = glht1$coefficients 
  estimates.j[2] = glht2$coefficients 
  estimates.j[3] = glht3$coefficients 
  estimates.j[4] = glht4$coefficients 
  p_values.j[1] = glht1$pvalues
  p_values.j[2] = glht2$pvalues
  p_values.j[3] = glht3$pvalues
  p_values.j[4] = glht4$pvalues
  estimates = rbind(estimates, estimates.j)
  p_values = rbind(p_values, p_values.j)
} 

dim(p_values)   
head(p_values)
alpha = 0.05 
sum(p_values[,1] < alpha)
sum(p_values[,2] < alpha)
sum(p_values[,3] < alpha) 
sum(p_values[,4] < alpha)  

p_value.1 = as.numeric(p.adjust(p_values[,1], "BH") < alpha)
p_value.2 = as.numeric(p.adjust(p_values[,2], "BH") < alpha)
p_value.3 = as.numeric(p.adjust(p_values[,3], "BH") < alpha)
p_value.4 = as.numeric(p.adjust(p_values[,4], "BH") < alpha)
sum(p_value.1) 
sum(p_value.2) 
sum(p_value.3)  
sum(p_value.4)   

sigmap.1 = sigmap.2 = sigmap.3 =sigmap.4 = mean.matrix.1 = mean.matrix.2 = mean.matrix.3 = mean.matrix.4 = matrix(NA, nrow=p,ncol=p) 
var.names = paste0("N",1:p) 
dimnames(sigmap.1) = dimnames(sigmap.2)=dimnames(sigmap.3)=dimnames(sigmap.4) =dimnames(mean.matrix.1) = dimnames(mean.matrix.2)=dimnames(mean.matrix.3)=dimnames(mean.matrix.4) = list(var.names,var.names) 
diag_indices = cbind(1:nrow(sigmap.1), 1:ncol(sigmap.1))

sigmap.1[lower.tri(sigmap.1)] = p_value.1[1:ncol(FZ.corr)]
sigmap.2[lower.tri(sigmap.2)] = p_value.2[1:ncol(FZ.corr)]
sigmap.3[lower.tri(sigmap.3)] = p_value.3[1:ncol(FZ.corr)] 
sigmap.4[lower.tri(sigmap.4)] = p_value.4[1:ncol(FZ.corr)]  
sigmap.1[diag_indices] = p_value.1[ncol(FZ.corr)+1:ncol(logvar.diag)]
sigmap.2[diag_indices] = p_value.2[ncol(FZ.corr)+1:ncol(logvar.diag)]
sigmap.3[diag_indices] = p_value.3[ncol(FZ.corr)+1:ncol(logvar.diag)]
sigmap.4[diag_indices] = p_value.4[ncol(FZ.corr)+1:ncol(logvar.diag)]

mean.matrix.1[lower.tri(sigmap.1)] = estimates[1:ncol(FZ.corr),1]
mean.matrix.2[lower.tri(sigmap.2)] = estimates[1:ncol(FZ.corr),2]
mean.matrix.3[lower.tri(sigmap.3)] = estimates[1:ncol(FZ.corr),3]
mean.matrix.4[lower.tri(sigmap.4)] = estimates[1:ncol(FZ.corr),4]
mean.matrix.1[diag_indices] = estimates[ncol(FZ.corr)+1:ncol(logvar.diag),1]
mean.matrix.2[diag_indices] = estimates[ncol(FZ.corr)+1:ncol(logvar.diag),2]
mean.matrix.3[diag_indices] = estimates[ncol(FZ.corr)+1:ncol(logvar.diag),3]
mean.matrix.4[diag_indices] = estimates[ncol(FZ.corr)+1:ncol(logvar.diag),4]


###################################################
#### Draw plots of the element-wise regression ####
###################################################
sigmap0 = sigmap.1  ## change it to sigmap.2, sigmap.3, sigmap.4
mean.matrix0 = mean.matrix.1  ## change it to mean.matrix.2, mean.matrix.3, mean.matrix.4

## a) significance matrix  
## Melt the correlation matrix   
melted_cormat = melt(sigmap0, na.rm = TRUE)
## Heatmap 
plot.sigmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(
    low = "blue", high = "green", mid = "white",  
    space = "Lab",  
    name=" ") + xlab(" ") + ylab(" ") + theme(
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      legend.key = element_rect(fill = "white")) + 
  guides(fill = guide_legend(override.aes= list(alpha = 0, color = "white"))) +
  theme(legend.key=element_rect(colour="white")) 
plot.sigmap 

## b) mean matrix  
## Get lower triangle of the correlation matrix
get_lower_tri=function(cormat){
  cormat[upper.tri(cormat)] = NA
  return(cormat)
}
lower_tri = get_lower_tri(mean.matrix0)

# Melt the correlation matrix 
melted_cormat = melt(lower_tri, na.rm = TRUE)
# Heatmap 
plot.mean = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",  
                       space = "Lab", 
                       #name="Pearson\nCorrelation") +
                       name=" ") + xlab(" ") + ylab(" ")
plot.mean 
 
grid.arrange(plot.sigmap, plot.mean, nrow=1, 
             top=textGrob("Covariance contrast, comparing 1) Short vs. Conventional sleeper, among Male", gp=gpar(fontsize=15)))
#grid.arrange(plot.sigmap, plot.mean, nrow=1, 
#             top=textGrob("Covariance contrast, comparing 2) Short vs. Conventional sleeper, among Female", gp=gpar(fontsize=15)))
#grid.arrange(plot.sigmap, plot.mean, nrow=1, 
#             top=textGrob("Covariance contrast, comparing 3) Male vs. Female, among Short sleeper", gp=gpar(fontsize=15)))
#grid.arrange(plot.sigmap, plot.mean, nrow=1, 
#             top=textGrob("Covariance contrast, comparing 4) Male vs. Female, among Conventional sleeper", gp=gpar(fontsize=15)))




 
######################
######################
### Now, fit BCAP ####
######################
######################
dat$nt = table(data$Subject)    
Y = do.call(rbind, dat$Yc) %>% as.matrix   
S_bar = dat$S_bar = apply(dat$S, c(2,3), mean) 
svd_S_bar = svd(S_bar)
whitening_mat = svd_S_bar$u %*% diag(1/sqrt(svd_S_bar$d)) %*% t(svd_S_bar$v) 
whitening_mat  
t_vec = sapply(dat$Yc, nrow)   
X.mat = model.matrix(as.formula("~sleep_duration + Gender + sleep_duration:Gender"), data=as.data.frame(do.call(rbind, dat$x)))
X = model.matrix(as.formula("~sleep_duration + Gender + sleep_duration:Gender"), data=data)[,-1]
attr(X, "assign") = NULL
colnames(Y)   

stan_data =list(p = ncol(Y),
                q = ncol(X),
                N = nrow(X),
                n = length(subject_ids),
                id= as.numeric(factor(data$Subject, levels=unique(data$Subject))), 
                X = X,
                Y = Y,  
                S= dat$S,
                whitening_mat=whitening_mat, 
                t_vec = t_vec,  # number of time points 
                Y.list = dat$Yc,
                X.mat= X.mat, 
                beta0_sd=2.5,
                beta_sd=2.5,
                eta=1)

## change the "iter_warmup" and "iter_sampling" below appropriately
## iter_warmup  = 700, iter_sampling= 1300
res = res_s1 = bcap_estimation(stan_data, init="cap", d=4, iter_warmup= 100, iter_sampling= 200)  
#load("res_03142024.RData")


#############################
## write a contrast vector ##
#############################
K1 = matrix(c(0, 1, 0, 1), 1) #  1) short vs. conventional sleeper among male
K2 = matrix(c(0, 1, 0, 0), 1) #  2) short vs. conventional sleeper among female;
K3 = matrix(c(0, 0, 1, 1), 1) #  3) male vs. female among short sleeper; and
K4 = matrix(c(0, 0, 1, 0), 1) #  4) male vs. female among conventional sleeper.

contrast.1 = compute_contrast(res, contr.vec = K1) #  1) short vs. conventional sleeper among male
contrast.2 = compute_contrast(res, contr.vec = K2) #  2) short vs. conventional sleeper among female; 
contrast.3 = compute_contrast(res, contr.vec = K3) #  3) male vs. female among short sleeper; and 
contrast.4 = compute_contrast(res, contr.vec = K4) #  4) male vs. female among conventional sleeper. 



################################
### 1) Variance Ratio plots  ###
################################
contrast = contrast.1  ## can change to contrast.2, contrast.3, contrast.4 
beta = contrast$log.cov.diag
#dim(beta)
#[1] 1300   15 
plot_brain = function (brain_slice, title, size=3) 
{
  img = readPNG(brain_slice)
  m = 10
  d = dim(img)
  img = img[m:(d[1] - m), m:(d[2] - m), ]
  d = dim(img)
  ggplot() + annotation_raster(img, xmin = -Inf, xmax = Inf, 
                               ymin = -Inf, ymax = Inf) + geom_point() + coord_fixed(ratio = d[1]/d[2]) + 
    annotate("text", x = 0, y = 0, label = title, size = size,#8, 
             color = "white") + theme_void()
}


### adopted from "CovRegFC::plot_coeff" with minor modifications 
plot_coeff_cap = function(beta, response, alpha = 0.05, col_number =1,
                          title = "", subtitle = "",brain_slices = NULL,  
                          coeff_labels = response,
                          lower.lim=0.5,upper.lim=2,
                          brain_nrow=3, 
                          thresh=0.02,
                          exponentiate=TRUE,
                          xintercept=1) 
{
  
  num_sim = dim(beta)[1]
  beta_2_perc = apply(beta, MARGIN = 2, 
                      function(vec) quantile(vec, probs = c(alpha/2, 0.5, 1 - alpha/2)))
  df_combo = data.frame(x=factor(response,levels=response), t(beta_2_perc))
  names(df_combo) = c("channel_name","low","mid","high")
  
  df_combo$confidence = abs(df_combo$mid)/abs(df_combo$high-df_combo$low)
  detect_sign = function(a, b, c) {
    if (sign(a) > 0 & sign(b) > 0 & abs(c)>thresh) {
      return("positive")
    }
    else if (sign(a) < 0 & sign(b) < 0 & abs(c)>thresh) {
      return("negative")
    }
    else {
      return("unclear")
    }
  }
  df_combo$sign = sapply(1:nrow(df_combo), 
                         function(i) detect_sign(df_combo$low[i],df_combo$high[i],df_combo$mid[i])) %>% factor(., levels = c("positive", "negative", "unclear"))
  df_combo$channel_name = factor(df_combo$channel_name, levels = df_combo$channel_name)
  beta_da_pos = lapply(1:num_sim, function(i) beta[i,] < 0) %>% Reduce("+", .)
  beta_da_neg = lapply(1:num_sim, function(i) beta[i,] > 0) %>% Reduce("+", .)
  ids = c(beta_da_neg, beta_da_pos)/num_sim <= alpha 
  neg_ids = (df_combo$sign=="negative") %>% which
  pos_ids = (df_combo$sign=="positive") %>% which
  
  df_combo = df_combo %>% mutate(color = ifelse(sign=="positive", "red", "blue"))
  df_combo$color[df_combo$sign=="unclear"] = "gray"
  
  if(exponentiate) df_combo[,c("low","mid","high")] = t(exp(beta_2_perc))
  if(is.null(lower.lim)) lower.lim = min(df_combo[,"low"])
  if(is.null(upper.lim)) upper.lim = min(df_combo[,"high"]) 
  
  if (!is.null(brain_slices)) {
    p_coeff = ggplot(df_combo, aes(x = mid, y = channel_name)) + 
      geom_vline(xintercept = xintercept, col = "black") + geom_point(size=1.5, colour=df_combo$color) + 
      geom_segment(mapping = aes(x = low, y = channel_name, 
                                 xend = high, yend = channel_name),
                   colour=df_combo$color) + labs(title = title, subtitle=subtitle) +  
      theme(axis.title.y = element_blank()) + 
      scale_y_discrete(breaks = coeff_labels)
    if(exponentiate) p_coeff  = p_coeff+ scale_x_continuous(trans = "log",limits=c(lower.lim, upper.lim))
    p_coeff
    
    ps = lapply(neg_ids, function(i) plot_brain(brain_slices[i],title = paste0("Net", i)))
    p_brains_neg = NULL
    if (length(ps) > 0){
      p_brains_neg = do.call(plot_grid, c(ps,ncol = 3,nrow = brain_nrow))  
    }
    
    ps = lapply(pos_ids, function(i) plot_brain(brain_slices[i],title = paste0("Net", i)))
    p_brains_pos = NULL
    if (length(ps) > 0){
      p_brains_pos = do.call(plot_grid, c(ps, ncol = 3,nrow = brain_nrow)) 
    } 
  }
  
  return(list(p_brains_neg=p_brains_neg, 
              p_coeff=p_coeff, 
              p_brains_pos=p_brains_pos,
              df_combo=df_combo))
}


response = paste("Net",sep="",1:p)
alpha = 0.05  
path = paste0("HCP_PTN820/groupICA/groupICA_3T_HCP820_MSMAll_d",num_regions,".ica/melodic_IC_sum.sum")
pngs = list.files(path = path,pattern = ".png")
brain_slices = full_paths = paste(path,pngs,sep = "/") 


contrast = contrast.1  ## can change to contrast.2, contrast.3, contrast.4 
p_brain_coeff_1 = plot_coeff_cap(beta=contrast$log.cov.diag, 
                                  response = response,
                                  alpha = 0.05,
                                  lower.lim=0.7,upper.lim=1.5, 
                                  title = "Variance Ratio (and 95% CrI), comparing",
                                  subtitle="1) Short vs. Conventional sleeper, among Male",
                                  brain_slices = full_paths) 
p_brain_coeff_1$p_brains_pos

contrast = contrast.2
p_brain_coeff_2 = plot_coeff_cap(beta=contrast$log.cov.diag, 
                                  response = channel_names,
                                  alpha = 0.05,
                                  lower.lim=0.7,upper.lim=1.5,
                                  title = "Variance Ratio (and 95% CrI), comparing",
                                  subtitle="2) Short vs. Conventional sleeper, among Female",
                                  brain_slices = full_paths) 

contrast = contrast.3
p_brain_coeff_3 = plot_coeff_cap(beta=contrast$log.cov.diag, 
                                  response = channel_names,
                                  alpha = 0.05,
                                  lower.lim=0.7,upper.lim=1.5,
                                  title = "Variance Ratio (and 95% CrI), comparing",
                                  subtitle="3) Male vs. Female, among Short sleeper",
                                  brain_slices = full_paths) 

contrast = contrast.4
p_brain_coeff_4 = plot_coeff_cap(beta=contrast$log.cov.diag, 
                                  response = channel_names,
                                  alpha = 0.05,
                                  lower.lim=0.7,upper.lim=1.5,
                                  title = "Variance Ratio (and 95% CrI), comparing",
                                  subtitle="4) Male vs. Female, among Conventional sleeper",
                                  brain_slices = full_paths) 


library(grid)
library(gridExtra) 
p_brain_coeff_1b = plot_grid(p_brain_coeff_1$p_coeff + xlab("Variance Ratio"),
                             p_brain_coeff_2$p_coeff + xlab("Variance Ratio"), 
                             p_brain_coeff_1$p_brains_pos, #+ggtitle("Parcel Set 2"),
                             labels = c(' ',
                                        ' ',
                                        'Short vs. Conventional Parcel Set'),
                             label_size = 10,
                             ncol = 3, align = "hv")
p_brain_coeff_1b
p_brain_coeff_2b = plot_grid(p_brain_coeff_3$p_coeff + xlab("Variance Ratio"),
                             p_brain_coeff_4$p_coeff + xlab("Variance Ratio"), 
                             p_brain_coeff_3$p_brains_pos, #+ggtitle("Parcel Set 2"),
                             labels = c(' ',
                                        ' ',
                                        'Male vs. Female Parcel Set'),
                             label_size = 10,
                             ncol = 3, align = "hv")
p_brain_coeff_2b  
grid.arrange(p_brain_coeff_1b, p_brain_coeff_2b, nrow=2)



 
#######################################################################
### 2) Covariance Contrast plots
#######################################################################

var.names = paste0("N",1:p) 
contrast = contrast.1

# a) significance matrix 
CrI_offdiag = t(apply(contrast$log.cov.offdiag, 2, 
                       function(x) quantile(x, probs =c(0.025,0.5,0.975))))
CrI_diag = t(apply(contrast$log.cov.diag, 2, 
                    function(x) quantile(x, probs =c(0.025,0.5,0.975))))
sigmap = matrix(NA, nrow=p,ncol=p)
diag_indices = cbind(1:nrow(sigmap), 1:ncol(sigmap)) 
sigmap[lower.tri(sigmap.1)] = as.numeric((sign(CrI_offdiag[,1])*sign(CrI_offdiag[,3])>0)&(abs(CrI_offdiag[,2])>0.02))  
sigmap[diag_indices] = as.numeric((sign(CrI_diag[,1])*sign(CrI_diag[,3])>0)&(abs(CrI_diag[,2])>0.02))
dimnames(sigmap) = list(var.names,var.names) 

## Melt the correlation matrix 
melted_cormat = melt(sigmap, na.rm = TRUE)
## Heatmap 
plot.sigmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",  
                       space = "Lab",  
                       name=" ") + xlab(" ") + ylab(" ") + theme(
                         legend.text = element_text(color = "white"),
                         legend.title = element_text(color = "white"),
                         legend.key = element_rect(fill = "white")) + 
  guides(fill = guide_legend(override.aes= list(alpha = 0, color = "white"))) +
  theme(legend.key=element_rect(colour="white")) 
plot.sigmap 

## b) mean matrix 
mean.matrix = apply(contrast$log.cov, c(2,3), mean)  
dimnames(mean.matrix) = list(var.names,var.names) 
## Get lower triangle of the correlation matrix
get_lower_tri=function(cormat){
  cormat[upper.tri(cormat)] = NA
  return(cormat)
}
lower_tri = get_lower_tri(mean.matrix)

## Melt the correlation matrix 
melted_cormat = melt(lower_tri, na.rm = TRUE)
## Heatmap 
plot.mean = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",  
                       space = "Lab",  
                       name=" ") + xlab(" ") + ylab(" ")
plot.mean 

grid.arrange(plot.sigmap, plot.mean, nrow=1, 
             top=textGrob("Covariance contrast, comparing 1) Short vs. Conventional sleeper, among Male", gp=gpar(fontsize=15)))
grid.arrange(plot.sigmap, plot.mean, nrow=1, 
             top=textGrob("Covariance contrast, comparing 2) Short vs. Conventional sleeper, among Female", gp=gpar(fontsize=15)))
grid.arrange(plot.sigmap, plot.mean, nrow=1, 
             top=textGrob("Covariance contrast, comparing 3) Male vs. Female, among Short sleeper", gp=gpar(fontsize=15)))
grid.arrange(plot.sigmap, plot.mean, nrow=1, 
             top=textGrob("Covariance contrast, comparing 4) Male vs. Female, among Conventional sleeper", gp=gpar(fontsize=15)))






#####################################################
### 3) display the estimated regression paramters ###
#####################################################
library(posterior)
library(bayesplot)  
ls(res) 
print(names(res$fit)) 
## a) Regression coefficients (B)
P= plot(res$fit,  pars =c("beta[2,1]", "beta[2,2]","beta[2,3]",  
                          "beta[1,1]", "beta[1,2]","beta[1,3]",  
                          "beta[3,1]", "beta[3,2]","beta[3,3]", 
                          "beta[4,1]", "beta[4,2]","beta[4,3]"), 
        point_est = "median",
        ci_level = 0.95, outer_level = 0.95, 
        fill_color = rep(c("maroon","maroon","maroon", "plum4","plum4","plum4"),2), 
        prob = 0.5, prob_outer = 0.95)

P + labs(
  title = "Regression coefficients (B)" 
) + vline_0() + scale_y_discrete(
    labels = c("beta[2,1]" = "C1 Sleep duration (short)",
               "beta[2,2]" = "C1 Gender (male)",
               "beta[2,3]" = "C1 Sleep x Gender Interaction",
               "beta[1,1]" = "C2 Sleep duration (short)",
               "beta[1,2]" = "C2 Gender (male)",
               "beta[1,3]" = "C2 Sleep x Gender Interaction",
               "beta[3,1]" = "C3 Sleep duration (short)",
               "beta[3,2]" = "C3 Gender (male)",
               "beta[3,3]" = "C3 Sleep x Gender Interaction",
               "beta[4,1]" = "C4 Sleep duration (short)",
               "beta[4,2]" = "C4 Gender (male)",
               "beta[4,3]" = "C4 Sleep x Gender Interaction"),
    limits = c("beta[4,3]", "beta[4,2]", "beta[4,1]",
               "beta[3,3]", "beta[3,2]", "beta[3,1]",
               "beta[1,3]", "beta[1,2]", "beta[1,1]",
               "beta[2,3]", "beta[2,2]", "beta[2,1]")
  ) + theme_gray(base_size = 16) + theme(axis.title.y=element_blank(),
                                         axis.title.x=element_blank()) #+ xlim(c(-0.05,0.55))


## b) random effect covariance parameter 
P2 = plot(res$fit,  pars =c("Omega[2,2]","Omega[1,1]", "Omega[3,3]", "Omega[4,4]",
                            "Omega[2,1]","Omega[2,3]","Omega[2,4]", 
                            "Omega[1,3]","Omega[1,4]","Omega[3,4]"), 
          point_est = "median",
          ci_level = 0.95, outer_level = 0.95, 
          fill_color = c(rep("royalblue1",10)), 
          prob = 0.5, prob_outer = 0.95)

P2 + labs(
  title = "Random effect covariance parameters" 
) + scale_y_discrete(
    labels = c("Omega[2,2]" = "Omega11",
               "Omega[1,1]" = "Omega22",
               "Omega[3,3]" = "Omega33",
               "Omega[4,4]" = "Omega44",
               "Omega[2,1]" = "Omega12",
               "Omega[2,3]" = "Omega13",
               "Omega[2,4]" = "Omega14",
               "Omega[1,3]" = "Omega23",
               "Omega[1,4]" = "Omega24",
               "Omega[3,4]" = "Omega34"),
    limits = c("Omega[3,4]", "Omega[1,4]", "Omega[1,3]",
               "Omega[2,4]", "Omega[2,3]", "Omega[2,1]",
               "Omega[4,4]", "Omega[3,3]", "Omega[1,1]",
               "Omega[2,2]")
  ) + theme_gray(base_size = 16) + theme(axis.title.y=element_blank(),
                                         axis.title.x=element_blank()) #+ xlim(c(0,1.12))



###########################################
### 4) display the loading coefficients ###
###########################################
contrast = res$draws_Gamma
dim(res$draws_Gamma) 
res$draws_Gamma[,,1]
p_gamma_coeff_1 = plot_coeff_cap(beta=res$draws_Gamma[,,1], 
                                  response = response,
                                  alpha = 0.05,#/p, 
                                  subtitle = "C1 loading coefficients", 
                                  brain_slices = full_paths,
                                  exponentiate=FALSE,
                                  xintercept=0)
p_gamma_coeff_1$p_coeff + xlab("coefficient") 

p_gamma_coeff_2 = plot_coeff_cap(beta=res$draws_Gamma[,,2], 
                                  response = response,
                                  alpha = 0.05,#/p, 
                                  subtitle = "C2 loading coefficients", 
                                  brain_slices = full_paths,
                                  exponentiate=FALSE,
                                  xintercept=0)
p_gamma_coeff_2$p_coeff + xlab("coefficient")


p_gamma_coeff_3 = plot_coeff_cap(beta=res$draws_Gamma[,,3], 
                                  response = response,
                                  alpha = 0.05,#/p, 
                                  subtitle = "C3 loading coefficients", 
                                  brain_slices = full_paths,
                                  exponentiate=FALSE,
                                  xintercept=0)
p_gamma_coeff_3$p_coeff + xlab("coefficient")


p_gamma_coeff_4 = plot_coeff_cap(beta=res$draws_Gamma[,,4], 
                                  response = response,
                                  alpha = 0.05,#/p, 
                                  subtitle = "C4 loading coefficients", 
                                  brain_slices = full_paths,
                                  exponentiate=FALSE,
                                  xintercept=0)  
p_gamma_coeff_4$p_coeff + xlab("coefficient") 

grid.arrange(p_gamma_coeff_1$p_coeff + xlab("coefficient"),
             p_gamma_coeff_2$p_coeff + xlab("coefficient"),
             p_gamma_coeff_3$p_coeff + xlab("coefficient"),
             p_gamma_coeff_4$p_coeff + xlab("coefficient"), 
             nrow=2)



###########################
### the end of the code ###
###########################