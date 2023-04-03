library(truncnorm)

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
    mle_model <- glm(Y[c,] ~ X[c,,] + Z[c,,] - 1, family = 'poisson', offset = rep(D[c], n_days))
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
  mle_model <- 
    glm(Y_pooled ~ X_pooled + Z_pooled - 1, family = 'poisson', offset = D_pooled)
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
    
    
  }
  
  
  ### Burn the Burn in 
  theta <- theta[-c(1:burn_in),,]
  mu <- mu[-c(1:burn_in),]
  eta <- eta[-c(1:burn_in),,]
  gamma <- gamma[-c(1:burn_in),]
  
  ### Package into list and return
  draws <- 
    list('theta' = theta,
         'mu' = mu,
         'eta' = eta,
         'gamma' = gamma)
  
  return(draws)
  
  
}