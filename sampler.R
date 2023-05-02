library(truncnorm)
library(mvtnorm)

### Compute Beta via profile likelihood approach
### This involves using theta and X in the offset to only estimate the 
### parameters for Z (beta)
update_beta <- function(theta, X, Y, Z, D) {
  beta <- matrix(NA, nrow = dim(X)[1], ncol = dim(Z)[3])
  for(c in 1:(dim(X)[1])) {
    mle_model <- 
      glm(Y[c,] ~ Z[c,,] - 1, 
          family = 'poisson', 
          offset = X[c,,] %*% theta[c,] + log(rep(D[c], (dim(X)[2]))))
    beta[c,] <- mle_model$coefficients
    
  }
  return(beta)
}

### Get the 2nd half of the coefficients from the profile liklihood 
### This involves using beta and Z in the offset to only estimate the 
### parameters for X (theta)
hat_theta <- function(beta, X, Y, Z, D) {
  theta_hat <- matrix(NA, nrow = dim(X)[1], ncol = dim(X)[3])
  Sigma <- list()
  
  for(c in 1:(dim(X)[1])) {
    mle_model <- 
      glm(Y[c,] ~ X[c,,] - 1, 
          family = 'poisson', 
          offset = Z[c,,] %*% beta[c,] + log(rep(D[c], (dim(X)[2]))))
    theta_hat[c,] <- mle_model$coefficients
    Sigma[[c]] <- vcov(mle_model)
  }
  
  return(list('theta_hat' = theta_hat,
              'Sigma' = Sigma))
}

### This proposes theta from Peng et. al rigorously derived proposal distribution
propose_theta <- function(mu, B, sigma2_gamma, Omega_gamma, theta_hat) {
  theta_prop <- matrix(NA, nrow = nrow(theta_hat), ncol = ncol(theta_hat))
  mean_vec <- map(1:nrow(theta_prop), ~{mu + B[[.x]] %*% (theta_hat[.x,]- mu)})
  var_mat <-  map(1:nrow(theta_prop), ~{sigma2_gamma * (diag(ncol(B[[.x]])) - B[[.x]]) %*% Omega_gamma })
  var_mat <- map(var_mat, ~{ 0.5 * (.x + t(.x)) + diag(seq(0.001, 0.0001, length.out = nrow(.x)), nrow = nrow(.x), ncol = ncol(.x))})
  
  for(c in 1:nrow(theta_prop)) {
    theta_prop[c,] <- rmvnorm(n = 1, mean = mean_vec[[c]], sigma = var_mat[[c]])
  }
  
  return(theta_prop)
}

### Accept or Reject Theta[c] samples
accept_reject_theta <- function(theta, theta_prop, X, Y, mu, Omega_gamma, sigma2_gamma) {
  alpha <- rep(NA, dim(X)[1])
  
  for(c in 1:(dim(X)[1])) {
    ### Compute acceptance prob
    prior <- dmvnorm(x = theta[c,], mean = mu, sigma = sigma2_gamma * Omega_gamma, log = T)
    prior_prop <- dmvnorm(x = theta_prop[c,], mean = mu, sigma = sigma2_gamma * Omega_gamma, log = T)
    log_lik <- sum(dpois(Y[c,], lambda =  exp(X[c,,] %*% theta[c,]), log = T))
    log_lik_prop <- 
      sum(dpois(Y[c,], lambda =  exp(X[c,,] %*% theta_prop[c,]), log = T))
    
    alpha[c] <- min(exp(prior_prop - prior + log_lik_prop - log_lik), 1)
    
    ### Replace w/ prob alpha
    if(runif(1) <= alpha[c]) {
      theta[c,] <- theta_prop[c,] 
    }
  }
  
  return(theta)
  
}

### MH step for eta
accept_reject_eta <- function(mu, eta, eta_prop, sigma2_eta, Omega_eta, Omega_eta_prop, theta) {
  lik <- 0
  lik_prop <- 0
  
  for(c in 1:nrow(theta)) {
    lik <- lik + dmvnorm(x = theta[c,], mean = mu, sigma = Omega_eta * sigma2_eta, log = T) 
    lik_prop <- lik_prop + dmvnorm(x = theta[c,], mean = mu, sigma = Omega_eta_prop * sigma2_eta, log = T) 
  }
  
  alpha <- min(exp(lik_prop - lik), 1)
  if(runif(1) <= alpha) {
    eta <- eta_prop
  }
  
  return(eta)
}

### MH step for gamma
accept_reject_gamma <- function(mu, gamma, gamma_prop, sigma2_gamma, Omega_gamma, Omega_gamma_prop) {
  
  lik <-  dmvnorm(x = mu, mean = rep(0, ncol(Omega_gamma)), sigma = Omega_gamma * sigma2_gamma, log = T) 
  lik_prop <-  dmvnorm(x = mu, mean = rep(0, ncol(Omega_gamma)), sigma = Omega_gamma_prop * sigma2_gamma, log = T) 
  
  
  alpha <- min(exp(lik_prop - lik), 1)
  if(runif(1) <= alpha) {
    gamma <- gamma_prop
  }
  
  return(gamma)
}



####################
### Main Sampler ###
####################
gibbs_sampler <- function(X, Y, Z, D, 
                          n_samples, 
                          burn_in, 
                          eta_hyperpriors = c(0.2, 0.8), 
                          gamma_hyperpriors = c(0.05, 0.75)) {
  
  ### Initialize Storage for parameters mu, theta, beta, gamma, eta, and sigmas 
  n_counties <- dim(X)[1]
  n_days <- dim(X)[2]
  n_lags <- dim(X)[3]
  n_iter <- n_samples + burn_in
  
  ### Beta
  ### 2-D Matrix where rows represent county and cols are the covariates we have 
  ### NOTE: We don't need to store Beta for each iteration because we don't 
  ### really care about doing any sort of inference for beta. We just need 
  ### it's current value
  beta <- 
    matrix(data = NA, 
           nrow = n_counties,
           ncol = dim(Z)[3], ### Extra 1 is for intercept 
           dimnames = dimnames(Z)[c(1,3)])
  
  
  ### Theta
  ### 3-D array where first dim is iteration, and next 2 are a matrix of n_counties x n_lags
  theta <- 
    array(data = NA, 
          dim = c(n_iter, n_counties, n_lags), 
          dimnames = list(1:n_iter, dimnames(X)[[1]], dimnames(X)[[3]]))
  
  ### Mu
  ### 2-D matrix where first dim (rows) is iteration and 2nd dim (cols) is n_lags
  mu <- 
    matrix(data = NA, 
           ncol = n_lags,
           nrow = n_iter,
           dimnames = list(1:n_iter, dimnames(X)[[3]]))
  
  ### Eta
  ### 2-D matrix where first dim (rows) is iteration and 2nd dim (cols) is eta_1 or eta_2
  eta <- 
    matrix(data = NA, 
           nrow = n_iter, 
           ncol = 2, 
           dimnames = list(1:n_iter, paste0('eta_', 1:2)))
  
  ### Sigma^2 eta
  ### (just 1 number that wil be fixed)
  sigma2_eta <- NA
  
  ### Gamma
  ### 2-D matrix where first dim (rows) is iteration and 2nd dim (cols) is gamma_1 or gamma_2
  gamma <- 
    matrix(data = NA, 
           nrow = n_iter, 
           ncol = 2, 
           dimnames = list(1:n_iter, paste0('gamma_', 1:2)))
  
  
  ### Sigma^2 gamma 
  ### (just 1 number that wil be fixed)
  sigma2_gamma <- NA
  
  
  ### Initialze Parameters
  for(c in 1:n_counties) {
    ### Beta and Theta
    ### Take MLE from Poisson GLM
    ### No intercept since I already included 1 column in the Z
    mle_model <- glm(Y[c,] ~ X[c,,] + Z[c,,] - 1, 
                     family = poisson(link = 'log'), 
                     offset = log(rep(D[c], n_days)))
    coeff <- mle_model$coefficients
    
    ### Theta
    theta[1,c, ] <- unname(coeff[1:n_lags])
    
    ### Beta
    beta[c,] <- unname(coeff[-c(1:n_lags)])
    
    ### Build X, Y, Z, pooled
    if(c == 1) {
      X_pooled <- X[c,,]
      Y_pooled <- Y[c,]
      Z_pooled <- Z[c,,]
      D_pooled <- rep(D[c], n_days)
    } else {
      X_pooled <- rbind(X_pooled, X[c,,])
      Y_pooled <- c(Y_pooled, Y[c,])
      Z_pooled <- rbind(Z_pooled, Z[c,,])
      D_pooled <- c(D_pooled, rep(D[c], n_days))
    }
  }
  
  ### Mu (MLE from pooled poisson regression)
  ### No convergence right now unless I get rid of offset
  mle_model <- 
    glm(Y_pooled ~ X_pooled + Z_pooled - 1, family = 'poisson', offset = log(D_pooled))
  mu[1,] <- mle_model$coefficients[1:n_lags]
  
  ### Gamma (Random Draw from Prior)
  gamma[1,] <- runif(2, min = gamma_hyperpriors[1], max = gamma_hyperpriors[2])
  
  ### Eta (Random Draw from Prior)
  eta[1,] <- runif(2, min = eta_hyperpriors[1], max = eta_hyperpriors[2])
  
  ### Sigma2 Gamma
  ### 10x the variance of the MLE of mu0
  sigma2_gamma <- 10 * summary(mle_model)$coefficients[1, 'Std. Error']^2
  
  ### Sigma2_eta
  ### 10x the variance of the MLE of mu0
  sigma2_eta <- 10 * summary(mle_model)$coefficients[1, 'Std. Error']^2
  
  
  ### Loop over iterations
  for(t in 2:n_iter) {
    ### Print every 100th iteration
    if(t %% 100 == 0) {
      cat('Iteration:', t, '\n')
    }
    ###################
    ### Update Beta ###
    ###################
    beta <- 
      update_beta(theta = theta[t-1,,],
                  X = X,
                  Y = Y,
                  Z = Z,
                  D = D)
    
    ####################
    ### Update Theta ###
    ####################
    ### Omega_gamma here is the covariance matrix (not scaled by sigma_eta)
    Omega_gamma <- 
      build_sigma(gamma1 = gamma[t-1, 1],
                  gamma2 = gamma[t-1, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    ### Compute theta hat c and corresponding variance matrix 
    tmp <- hat_theta(beta, X, Y, Z, D)
    theta_hat <- tmp$theta_hat
    Sigma <- tmp$Sigma
    
    ### Compute B_c matrix for each C
    B <- map(Sigma, ~{(sigma2_gamma *  Omega_gamma) * solve(.x + (sigma2_gamma *  Omega_gamma))})
    
    ### Propose Theta
    theta_prop <- 
      propose_theta(mu = mu[t-1,], 
                    B = B, 
                    sigma2_gamma = sigma2_gamma, 
                    Omega_gamma = Omega_gamma,
                    theta_hat = theta_hat)
    
    ### Decide to accept or reject each theta_c
    theta[t,,] <- 
      accept_reject_theta(theta = theta[t-1,,],
                          theta_prop = theta_prop,
                          X = X, 
                          Y = Y, 
                          mu = mu[t-1,],
                          Omega_gamma = Omega_gamma,
                          sigma2_gamma = sigma2_gamma)
    
    ####################
    ### Update Eta  ###
    ####################
    Omega_eta <- 
      build_sigma(gamma1 = eta[t-1, 1],
                  gamma2 = eta[t-1, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    ### Propose Eta 
    eta1_prop <- runif(1, min = max(eta_hyperpriors[1], eta[t-1, 1] - 0.005), max = min(eta_hyperpriors[2], eta[t-1, 1] + 0.005))
    eta2_prop <- runif(1, min = max(eta_hyperpriors[1], eta[t-1, 2] - 0.005), max = min(eta_hyperpriors[2], eta[t-1, 2] + 0.005))
    
    ### Accept/Reject Eta 1/2 (seperately)
    ### Eta 1
    Omega_prop1 <- 
      build_sigma(gamma1 =  eta1_prop,
                  gamma2 = eta[t-1, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    eta[t,] <- 
      accept_reject_eta(mu = mu[t-1,],
                        eta = eta[t-1,], 
                        eta_prop = c(eta1_prop, eta[t-1,2]),
                        sigma2_eta = sigma2_eta,
                        Omega_eta = Omega_eta, 
                        Omega_eta_prop =  Omega_prop1,
                        theta = theta[t,,])
    
    ### Eta 2
    Omega_prop2 <- 
      build_sigma(gamma1 =  eta[t, 2],
                  gamma2 = eta2_prop,
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    eta[t,] <- 
      accept_reject_eta(mu = mu[t-1,],
                        eta = eta[t,], 
                        eta_prop = c(eta[t,1], eta2_prop),
                        sigma2_eta = sigma2_eta,
                        Omega_eta = Omega_eta, 
                        Omega_eta_prop =  Omega_prop2,
                        theta = theta[t,,])
    
    
    ####################
    ### Update Gamma ###
    ####################
    ### Propose Gamma
    gamma1_prop <- runif(1, min = max(gamma_hyperpriors[1], gamma[t-1, 1] - 0.005), max = min(gamma_hyperpriors[2], gamma[t-1, 1] + 0.005))
    gamma2_prop <- runif(1, min = max(gamma_hyperpriors[1], gamma[t-1, 2] - 0.005), max = min(gamma_hyperpriors[2], gamma[t-1, 2] + 0.005))
    
    ### Accept/Reject Gamma 1/2 (seperately)
    ### Gamma1
    Omega_prop1 <- 
      build_sigma(gamma1 =  gamma1_prop,
                  gamma2 = gamma[t-1, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    gamma[t,] <- 
      accept_reject_gamma(mu = mu[t-1,],
                          gamma = gamma[t-1,], 
                          gamma_prop = c(gamma1_prop, gamma[t-1,2]),
                          sigma2_gamma = sigma2_gamma,
                          Omega_gamma = Omega_gamma, 
                          Omega_gamma_prop =  Omega_prop1)
    
    ### Gamma 2
    Omega_prop2 <- 
      build_sigma(gamma1 =  gamma[t, 2],
                  gamma2 = gamma2_prop,
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    gamma[t,] <- 
      accept_reject_gamma(mu = mu[t-1,],
                          gamma = gamma[t,], 
                          gamma_prop = c(gamma[t,1], gamma2_prop),
                          sigma2_gamma = sigma2_gamma,
                          Omega_gamma = Omega_gamma, 
                          Omega_gamma_prop =  Omega_prop2)
    
    ####################
    #### Update Mu #####
    ####################
    ### Re-compute (yet again!!!) using most up to date values
    Omega_gamma <- 
      build_sigma(gamma1 = gamma[t, 1],
                  gamma2 = gamma[t, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    Omega_eta <- 
      build_sigma(gamma1 = eta[t, 1],
                  gamma2 = eta[t, 2],
                  n_lags = n_lags-1,
                  sigma2 = 1)
    
    M <- (sigma2_gamma * Omega_gamma) %*% solve((sigma2_gamma * Omega_gamma) + 1/n_counties * sigma2_eta * Omega_eta)
    theta_bar <- colMeans(theta[t,,])
    cov_mat <- (diag(rep(1, ncol(M))) - M)  %*% (sigma2_gamma * Omega_gamma)
    cov_mat <- 0.5 * (cov_mat + t(cov_mat))
    mu[t,] <- rmvnorm(n = 1, mean = M %*% theta_bar, sigma = cov_mat)
    
  }
  
  
  ### Burn the Burn in 
  theta <- theta[-c(1:burn_in),,]
  mu <- mu[-c(1:burn_in),]
  eta <- eta[-c(1:burn_in),]
  gamma <- gamma[-c(1:burn_in),]
  
  ### Package into list and return
  draws <- 
    list('theta' = theta,
         'mu' = mu,
         'eta' = eta,
         'gamma' = gamma)
  
  return(draws)
  
  
}
