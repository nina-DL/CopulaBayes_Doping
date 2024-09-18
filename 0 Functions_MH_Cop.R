# Load / install necessary packages ---------------------------------------
pkgs = c("bayestestR", "truncnorm")
missing.packages = pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
sapply(pkgs, require, character.only = TRUE)

# Clayton copula ----------------------------------------------------------

# Log-density of a Clayton
ldClay = function(u,v,theta){
  if(theta < 0){stop("Error: The parameter cannot be smaller than 0")}
  if((any(u>1)|any(u<0)|any(v>1)|any(v<0))){stop("Pseudo-obs do not have a valid range: should be between 0 and 1")}
  
  num = (theta+1)*(u*v)^theta
  den = (u^theta + v^theta - (u*v)^theta)^(1/theta+2)
  ld = log(num) - log(den)
  return(ld)
}

# Log-density of a Clayton
ldClay2 = function(u,v,theta){
  if(theta < 0){stop("Error: The parameter cannot be smaller than 0")}   
  if((any(u>1)|any(u<0)|any(v>1)|any(v<0))){stop("Pseudo-obs do not have a valid range: should be between 0 and 1")}
  
  s1 = log(theta+1)-(1+theta)*log(u*v)
  s2 = -(1/theta+2)*log(u^(-theta) + v^(-theta) - 1)
  ld = s1+s2
  return(ld)
}

# Log-likelihood of a Clayton
lLik_Clay = function(u,v,theta){
  if(theta < 0){stop("Clayton parameter cannot be smaller than 0")}   
  if((any(u>1)|any(u<0)|any(v>1)|any(v<0))){stop("Pseudo-obs do not have a valid range: should be between 0 and 1")}
  
  f_Clay = log(theta+1)-(1+theta)*log(u*v) - (1/theta+2)*log(u^(-theta) + v^(-theta) - 1)
  lLik = sum(f_Clay)
  return(lLik)
}

# Log-posterior of a Clayton
lpost_Clay = function(u, v, theta, prior_spec = NULL, par_prior1, par_prior2){
  
  if(is.null(prior_spec)){
    warning("Prior was not specified: default 'unif' prior is used")   
    prior_spec = "unif"
    d_prior = dunif(theta, par_prior1, par_prior2, log = T)
  }
  
  if(prior_spec == "truncnorm"){
    d_prior = log(dtruncnorm(theta, a=0, b=100, mean = par_prior1, sd = par_prior2))
  } else if(prior_spec == "lognorm"){
    d_prior = dlnorm(theta, meanlog = log(par_prior1^2 / sqrt(par_prior2^2 + par_prior1^2)), 
                     sdlog = sqrt(log(1 + (par_prior2^2 / par_prior1^2))), log = T)
  } else if(prior_spec == "unif"){
    d_prior = dunif(theta, par_prior1, par_prior2, log = T)
  }
  
  lpost = lLik_Clay(u, v, theta) + d_prior
  return(lpost)
}

#lpost_Clay(u, v, theta = 1, prior = "lognorm", par_prior1 = par_prior1, par_prior2 = par_prior2)

RWMH_Clay = function(M, u, v, surv = T, theta0, par_prop = .6,
                     prior_spec = "unif", par_prior1 = 0, par_prior2 = 50, burnin = 0.1){
  
  # This function implements a random-walk Metropolis Hastings to simulate posterior draws from the Clayton copula 
  # It assumes a uniform proposal, but prior can be chosen: unif, lognorm, truncnorm
  # ARGUMENTS OF THE FUNCTION:
  # M :: length of the chain
  # u,v :: marginals (pseudo-obs)
  # surv :: whether we should refer to a Clayton (surv = F) or survival Clayton (surv = T)
  # theta0 :: initial value of the chain (by default MLE)
  # par_prop :: parameter (variability) of the proposal distribution; can be used to calibrate the acceptance rate [25-30% is ok]
  # prior_spec :: prior distribution; can be one of the following: "unif" (default), "lognorm", "truncnorm"
  # par_prior1 :: parameter 1 (min / mean) of the prior distribution
  # par_prior2 :: parameter 2 (max / sd) of the prior distribution
  # burnin :: set the burnin proportion in [0,1]
  
  if(surv==T){
    u=1-u
    v=1-v
  }
  
  b = theta0
  acra = 0
  theta_draws = c()
  
  Mb = M + M*burnin
  for(m in 1:Mb){
  
    r = runif(1, -par_prop, par_prop)
    a = r + b
    
    num = lpost_Clay(u = u, v = v, theta = a, prior = prior_spec, par_prior1 = par_prior1, par_prior2 = par_prior2)
    den = lpost_Clay(u = u, v = v, theta = b, prior = prior_spec, par_prior1 = par_prior1, par_prior2 = par_prior2)
    
    # with a lognorm proposal:
    # a = rlnorm(1, b, par_prop)
    # num.prop = dlnorm(b, a, par_prop, log = T)
    # den.prop = dlnorm(a, b, par_prop, log = T)
    # Post_target = exp(num - den + num.prop - den.prop)
    
    Post_target = exp(num - den)
    
    if(Post_target>runif(1)){
      b = a
      acra = acra+1
    }
    theta_draws[m] = b
  }
  AR = acra/Mb
  res = list(theta_draws = theta_draws[(M*burnin+1):Mb], MAP = map_estimate(theta_draws[(M*burnin+1):Mb], na.rm = T), AR = AR)
  return(res)
}

# post_theta = RWMH_Clay(M = 10000, burnin = 0.1, u = U1, v = U2, surv = T, theta0 = est_copula$par, par_prop = .1, prior_spec = "lognorm", par_prior1 = est_copula$par, par_prior2 = 100)
# 
# hist(post_theta$theta_draws)
# mean(post_theta$theta_draws)
# abline(v = mean(post_theta$theta_draws), col = "red")
# post_theta$AR



# Frank copula ------------------------------------------------------------

ldfr<-function(u,v,a){
  fa<-log(a-a*exp(-a))-a*(u+v)-log((exp(-a)-1+ (exp(-a*u)-1)*(exp(-a*v)-1))^2)
  return(fa)
}

RWMH_Frank = function(M, u, v, theta0, par_prop = .6,
                      par_prior1 = 5, par_prior2 = 1, burnin = 0.1){
  
  # This function implements a random-walk Metropolis Hastings to simulate posterior draws from the Frank copula 
  # It assumes a uniform proposal, and a Gamma prior
  # ARGUMENTS OF THE FUNCTION:
  # M :: length of the chain
  # u,v :: marginals (pseudo-obs)
  # theta0 :: initial value of the chain (by default MLE)
  # par_prop :: parameter (variability) of the proposal distribution; can be used to calibrate the acceptance rate [25-30% is ok]
  # par_prior1 :: parameter 1 (shape) of the prior distribution
  # par_prior2 :: parameter 2 (rate) of the prior distribution
  # burnin :: set the burnin proportion in [0,1]
  
  b = theta0
  acra = 0
  theta_draws = c()
  
  Mb = M + M*burnin
  for(m in 1:Mb){
    
    r = runif(1, -par_prop, par_prop)
    a = r + b
    
    num = ldfr(u = u, v = v, a)
    den = ldfr(u = u, v = v, b)
    LH = exp(sum(num - den))
    G<-((a/b)^par_prior1)*exp(par_prior2*(b-a)) # prior gamma
    if(is.na(LH*G) == T){LH = 0; G = 0}
    if(LH*G>runif(1)){
      b = a
      acra = acra+1
    }
    theta_draws[m] = b
  }
  AR = acra/Mb
  res = list(theta_draws = theta_draws[(M*burnin+1):Mb], MAP = map_estimate(theta_draws[(M*burnin+1):Mb], na.rm = T), AR = AR)
  return(res)
}
