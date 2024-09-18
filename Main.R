# Load / install necessary packages ---------------------------------------
library(parallel)
pkgs = c("MASS", "copula", "VineCopula", "KScorrect", "dplyr", "RcppEigen", "bayestestR", "brms", "tidybayes", "rstan", "mclust", "mvtnorm")
missing.packages = pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
sapply(pkgs, require, character.only = TRUE)

library(bmixture)
library(sn)
library(ggplot2)
require(gridExtra)

# Hypersetup: load data (S1) and necessary function -----------------------

load("S1_data_500.RData")
# mod_S1 records the true data generation model
mod_S1

# Dimension of the dataset
myN = dim(S1_draws_ext[[1]])[1]

# load the MH function for RWMH - Clay
source("0 Functions_MH_Cop.R")


# Predictive --------------

# function to run posterior predictive inference with the univariate model
get_pred_draws_uv = function(Mpred = 1e3, post_samples, seed = 123){
  
  Mpost = length(post_samples[[1]])
  
  Y_t1_uv = matrix(NA, ncol=2, nrow = Mpred)
  set.seed(seed)
  for(i in 1:Mpred){
    
    idx = sample(c(1:Mpost), 1)
    
    # SOTA: univariate
    Hgb_mean = post_samples$mu[idx]
    Hgb_sd = post_samples$s[idx]
    OFFs_mean = c(post_samples$mu1[idx], post_samples$mu2[idx], post_samples$mu3[idx])
    OFFs_sd = c(post_samples$s1[idx], post_samples$s2[idx], post_samples$s3[idx])
    OFFs_w = c(post_samples$t1[idx], post_samples$t2[idx], post_samples$t3[idx])
    
    Y_t1_uv[i,1] = rnorm(1, mean = Hgb_mean, sd = Hgb_sd)
    Y_t1_uv[i,2] = KScorrect::rmixnorm(n=2, mean=OFFs_mean, sd=OFFs_sd, pro = OFFs_w)[1]
  }
  return(Y_t1_uv = Y_t1_uv)
}

# function to run posterior predictive inference with the multivariate Gaussian model
get_pred_draws_MVN = function(Mpred = 1e3, post_samples, seed = 123){
  
  Mpost = length(post_samples[[1]])
  
  Y_t1_mv = matrix(NA, ncol=2, nrow = Mpred)
  set.seed(seed)
  for(i in 1:Mpred){
    
    idx = sample(1:Mpost, 1)
    
    # Greci: multivariate Normal
    Hgb_mean = post_samples$mu1[idx]
    Hgb_var = post_samples$s1[idx]^2
    
    OFFs_mean = post_samples$mu2[idx]
    OFFs_var = post_samples$s2[idx]^2
    
    HgbOFFs_cov = post_samples$rho[idx]*sqrt(Hgb_var)*sqrt(OFFs_var)
    
    Sigma = matrix(c(Hgb_var, rep(HgbOFFs_cov, 2), OFFs_var), nrow = 2)
    
    Y_t1_mv[i,] = rmvnorm(1, mean = c(Hgb_mean, OFFs_mean), sigma = Sigma, checkSymmetry = F, method = "chol")
    
  }
  return(Y_t1_mv = Y_t1_mv)
}

# function to run posterior predictive inference with the copula-based model
get_pred_draws_cop = function(Mpred = 1e3, post_samples, seed = 123){
  
  Mpost = length(post_samples[[1]])
  
  Y_t1_Clay = matrix(NA, ncol=2, nrow = Mpred)
  set.seed(seed)
  for(i in 1:Mpred){
    
    idx = sample(c(1:Mpost), 1)
    
    Hgb_mean = post_samples$mu[idx]
    Hgb_sd = post_samples$s[idx]
    OFFs_mean = c(post_samples$mu1[idx], post_samples$mu2[idx], post_samples$mu3[idx])
    OFFs_sd = c(post_samples$s1[idx], post_samples$s2[idx], post_samples$s3[idx])
    OFFs_w = c(post_samples$t1[idx], post_samples$t2[idx], post_samples$t3[idx])
    
    # PROPOSED: MVT with correct copula
    Clay_Copula = claytonCopula(param = post_samples$thetaC[idx], dim = 2)
    Clay_Copula = rotCopula(Clay_Copula)
    
    Biv_distC = mvdc(copula = Clay_Copula,
                     margins = c("norm", "mixnorm"),
                     paramMargins = list(list(mean = Hgb_mean, sd = Hgb_sd),
                                         list(mean=OFFs_mean, sd=OFFs_sd, pro=OFFs_w)))
    
    Y_t1_Clay[i,] = rMvdc(1, Biv_distC)
  }
  return(Y_t1_Clay = Y_t1_Clay)
}


# Hypersetup for stan -----------------------

rstan_options(auto_write = TRUE)

Hgb_m1 = formula('Hgb ~ 1')
OFFs_m1 = formula('OFFs ~ 1')

mix <- mixture(gaussian, gaussian, gaussian)
prior_OFFs <- c(
  prior(normal(80, 10), Intercept, dpar = mu1),
  prior(normal(90, 10), Intercept, dpar = mu2),
  prior(normal(100, 10), Intercept, dpar = mu3)
)

mvt_m <- 
  bf(mvbind(Hgb, OFFs) ~ 1) +
  set_rescor(TRUE)


# Univariate Model -----------------------

Univ = function(idx = 1, seed = 1234, n_mcmc = 5000){
  S1_draws = S1_draws_ext[[idx]]
  set.seed(seed)
  
  #marg1
  Post_draws_Hgb_m1 = list()
  Hgb_m1_fit = brm(formula = Hgb_m1, data = S1_draws, family = "normal", iter = n_mcmc, 
                   chains = 2, init  = "random", cores = 1, thin = 5, seed = 123)
  Hgb_m1_fit = Hgb_m1_fit %>% spread_draws(b_Intercept, sigma)
  Post_draws_Hgb_m1$mu = Hgb_m1_fit$b_Intercept
  Post_draws_Hgb_m1$s = Hgb_m1_fit$sigma
  
  #marg2
  Post_draws_OFFs_m1 = list()
  OFFs_m1_fit = brm(formula = OFFs_m1, data = S1_draws, family = mix, prior = prior_OFFs,
                    iter = n_mcmc, chains = 2, cores = 1, init  = "random", thin = 5, seed = 123)
  OFFs_m1_fit = OFFs_m1_fit %>% spread_draws(b_mu1_Intercept, b_mu2_Intercept, b_mu3_Intercept,
                                             sigma1, sigma2, sigma3, theta1, theta2, theta3)
  Post_draws_OFFs_m1$mu1 = OFFs_m1_fit$b_mu1_Intercept
  Post_draws_OFFs_m1$mu2 = OFFs_m1_fit$b_mu2_Intercept
  Post_draws_OFFs_m1$mu3 = OFFs_m1_fit$b_mu3_Intercept
  Post_draws_OFFs_m1$s1 = OFFs_m1_fit$sigma1
  Post_draws_OFFs_m1$s2 = OFFs_m1_fit$sigma2
  Post_draws_OFFs_m1$s3 = OFFs_m1_fit$sigma3
  Post_draws_OFFs_m1$t1 = OFFs_m1_fit$theta1
  Post_draws_OFFs_m1$t2 = OFFs_m1_fit$theta2
  Post_draws_OFFs_m1$t3 = OFFs_m1_fit$theta3
  
  post_samples_uni = data.frame(mu = Post_draws_Hgb_m1$mu,
                                s = Post_draws_Hgb_m1$s,
                                mu1 = Post_draws_OFFs_m1$mu1,
                                mu2 = Post_draws_OFFs_m1$mu2,
                                mu3 = Post_draws_OFFs_m1$mu3,
                                s1 = Post_draws_OFFs_m1$s1,
                                s2 = Post_draws_OFFs_m1$s2,
                                s3 = Post_draws_OFFs_m1$s3,
                                t1 = Post_draws_OFFs_m1$t1,
                                t2 = Post_draws_OFFs_m1$t2,
                                t3 = Post_draws_OFFs_m1$t3)
  
  Res_pred = get_pred_draws_uv(Mpred = 1e3, post_samples = post_samples_uni)
  
  return(Res_pred)
}


# Multivariate Model with Copula -----------------------

Multiv_cop = function(idx = 1, seed = 1234, n_mcmc = 5000){
  S1_draws = S1_draws_ext[[idx]]
  set.seed(seed)
  
  #marg1
  Post_draws_Hgb_m1 = list()
  Hgb_m1_fit = brm(formula = Hgb_m1, data = S1_draws, family = "normal", iter = n_mcmc, 
                   chains = 2, init  = "random", cores = 1, thin = 5, seed = 123)
  Hgb_m1_fit = Hgb_m1_fit %>% spread_draws(b_Intercept, sigma)
  Post_draws_Hgb_m1$mu = Hgb_m1_fit$b_Intercept
  Post_draws_Hgb_m1$s = Hgb_m1_fit$sigma
  
  #marg2
  Post_draws_OFFs_m1 = list()
  OFFs_m1_fit = brm(formula = OFFs_m1, data = S1_draws, family = mix, prior = prior_OFFs,
                    iter = n_mcmc, chains = 2, cores = 1, init  = "random", thin = 5, seed = 123)
  OFFs_m1_fit = OFFs_m1_fit %>% spread_draws(b_mu1_Intercept, b_mu2_Intercept, b_mu3_Intercept,
                                             sigma1, sigma2, sigma3, theta1, theta2, theta3)
  Post_draws_OFFs_m1$mu1 = OFFs_m1_fit$b_mu1_Intercept
  Post_draws_OFFs_m1$mu2 = OFFs_m1_fit$b_mu2_Intercept
  Post_draws_OFFs_m1$mu3 = OFFs_m1_fit$b_mu3_Intercept
  Post_draws_OFFs_m1$s1 = OFFs_m1_fit$sigma1
  Post_draws_OFFs_m1$s2 = OFFs_m1_fit$sigma2
  Post_draws_OFFs_m1$s3 = OFFs_m1_fit$sigma3
  Post_draws_OFFs_m1$t1 = OFFs_m1_fit$theta1
  Post_draws_OFFs_m1$t2 = OFFs_m1_fit$theta2
  Post_draws_OFFs_m1$t3 = OFFs_m1_fit$theta3
  
  #copula
  Hgb_mle <- fitdistr(S1_draws[,1],"normal")
  OFFs_mle = mclust::Mclust(S1_draws[,2], G = 3)$param #modelNames = mclustModelNames("V"), 
  
  U1 = pnorm(S1_draws[,1], mean = Hgb_mle$estimate[1], sd = Hgb_mle$estimate[2])
  U2 = KScorrect::pmixnorm(S1_draws[,2], mean=OFFs_mle$mean, sd=sqrt(OFFs_mle$variance$sigmasq), pro=OFFs_mle$pro)
  
  copula_mle = BiCopSelect(U1, U2, familyset = 0:6)
  Clay_mle = claytonCopula(par = copula_mle$par, dim = 2)
  surv_Clay_mle = rotCopula(Clay_mle)
  
  # Get posterior draws from the copula model
  Mpost = length(Post_draws_Hgb_m1$mu)
  Mpost_cop = 100
  # Version 2: Copula est conditioned on posterior draws (uncertainty propagation: full) -- similar fashion to https://rivista-statistica.unibo.it/article/download/5827/5541/16987
  Clay_post_chains = matrix(NA, ncol = Mpost, nrow = Mpost_cop)
  for(i in 1:Mpost){
    Hgb_mean = Post_draws_Hgb_m1$mu[i]
    Hgb_sd = Post_draws_Hgb_m1$s[i]
    OFFs_mean = c(Post_draws_OFFs_m1$mu1[i], Post_draws_OFFs_m1$mu2[i], Post_draws_OFFs_m1$mu3[i])
    OFFs_sd = c(Post_draws_OFFs_m1$s1[i], Post_draws_OFFs_m1$s2[i], Post_draws_OFFs_m1$s3[i])
    OFFs_w = c(Post_draws_OFFs_m1$t1[i], Post_draws_OFFs_m1$t2[i], Post_draws_OFFs_m1$t3[i])
    
    U1_post = pnorm(S1_draws[,1], mean = Hgb_mean, sd = Hgb_sd)
    U2_post = KScorrect::pmixnorm(S1_draws[,2], mean=OFFs_mean, sd=OFFs_sd, pro=OFFs_w)
    
    MH_res = RWMH_Clay(M = Mpost_cop, u = U1_post, v = U2_post, theta0 = copula_mle$par)
    Clay_post_chains[,i] = MH_res$theta_draws
  }
  
  Post_draws_Clay = Clay_post_chains[,Mpost_cop]
  
  post_samples_cop = data.frame(mu = Post_draws_Hgb_m1$mu,
                                s = Post_draws_Hgb_m1$s,
                                mu1 = Post_draws_OFFs_m1$mu1,
                                mu2 = Post_draws_OFFs_m1$mu2,
                                mu3 = Post_draws_OFFs_m1$mu3,
                                s1 = Post_draws_OFFs_m1$s1,
                                s2 = Post_draws_OFFs_m1$s2,
                                s3 = Post_draws_OFFs_m1$s3,
                                t1 = Post_draws_OFFs_m1$t1,
                                t2 = Post_draws_OFFs_m1$t2,
                                t3 = Post_draws_OFFs_m1$t3,
                                thetaC = Post_draws_Clay)
  
  Res_pred = get_pred_draws_cop(Mpred = 1e3, post_samples = post_samples_cop)
}


# Multivariate Normal Model -----------------------

Multiv_N = function(idx = 1, seed = 1234, n_mcmc = 5000){
  S1_draws = S1_draws_ext[[idx]]
  set.seed(seed)
  
  Post_draws_mvt_m1 = list()
  mvt_fit = brm(mvt_m, data = S1_draws, iter = n_mcmc, chains = 2,
                init  = "random", cores = 1, thin = 5, seed = 123)
  mvt_fit = mvt_fit %>% spread_draws(b_Hgb_Intercept, b_OFFs_Intercept,
                                     sigma_Hgb, sigma_OFFs, rescor__Hgb__OFFs)
  Post_draws_mvt_m1$mu1 = mvt_fit$b_Hgb_Intercept
  Post_draws_mvt_m1$mu2 = mvt_fit$b_OFFs_Intercept
  Post_draws_mvt_m1$s1 = mvt_fit$sigma_Hgb
  Post_draws_mvt_m1$s2 = mvt_fit$sigma_OFFs
  Post_draws_mvt_m1$rho = mvt_fit$rescor__Hgb__OFFs
  
  post_samples_mvn = data.frame(mu1 = Post_draws_mvt_m1$mu1,
                                mu2 = Post_draws_mvt_m1$mu2,
                                s1 = Post_draws_mvt_m1$s1,
                                s2 = Post_draws_mvt_m1$s2,
                                rho = Post_draws_mvt_m1$rho)
  
  Res_pred = get_pred_draws_MVN(Mpred = 1e3, post_samples = post_samples_mvn)
}



# Running times -------------------------------------------------------

Sys.info()

devtools::install_github("olafmersmann/microbenchmarkCore")
devtools::install_github("olafmersmann/microbenchmark")

library(microbenchmark)

set.seed(123)
runn_J <- microbenchmark("Univariate Model" = Univ(),
                         "Multivariate Normal Model" = Multiv_N(),
                         "Multivariate Copula-based Model" = Multiv_cop(),
                         times = 50L,
                         control = list(warmup = 2)
)


boxplot(runn_J, log = F, xlab = "Model", ylab = "Time [seconds]") 


# Compute the HPR ---------------------------------------------------------------------

install.packages("devtools")
devtools::install_github("nina-DL/HDR2D")
library(HDR2D)

# get some posterior predictive draws
Res_pred = Univ()

# Estimate a 95% HDR with Nonparametric Copula
HDR.2d(sn = Res_pred, measure = "DE.NPCop", coverage_prob = 0.95, build_plot = T, cex = 0.7, pch = 10)

# evaluate running times for HPR

set.seed(123)
runn_H <- microbenchmark("HPR" = HDR.2d(sn = Res_pred, measure = "DE.NPCop", coverage_prob = 0.95, build_plot = F, cex = 0.7, pch = 10),
                         times = 50L, unit = "seconds",
                         control = list(warmup = 2)
)

boxplot(runn_H, log = F, xlab = "Model", ylab = "Time [seconds]") 


