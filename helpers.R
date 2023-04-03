### Function to augment a dataset with lags of a particular variable
create_lag_df <- function(df_fips, n_lags, lag_var = 'PM25_mean_pred') {
  ### Create the lag columns by mapping over the function lag for each lag period
  ### map_dfc binds resulting columns together w/ bind_cols()
  df_lags <- 
    map_dfc(0:n_lags, ~{
      df <- tibble('lag' = lag(df_fips[[ lag_var ]], n = .x))
      names(df) <- paste0('lag_', .x)
      df
    })
  
  ### Augment our df_fips w/ the new lags
  df_fips <- 
    df_fips %>% 
    bind_cols(df_lags)
  
  return(df_fips)
}


### Function to build covariance matrix 
### gamma1, gamma2 = parameters that control the decay of the lag function 
### sigma2 = variance of coeff for lag_0
### n_lags = n_lags + 1 = dim of matrix 
build_sigma <- function(gamma1, gamma2, sigma2, n_lags) {
  Sigma <- matrix(0, nrow = n_lags+1, ncol = n_lags+1)
  colnames(Sigma) <- paste0('lag_', 0:n_lags)
  rownames(Sigma) <- paste0('lag_', 0:n_lags)
  
  for(l1 in 0:n_lags) {
    for(l2 in 0:n_lags) {
      ### Diagonal Entries
      if(l1 == l2) {
        Sigma[l1+1, l2+1] <- sigma2 * exp(-gamma1 * l1)
      } else if(abs(l1 - l2) >= 1) { ### Off diagonals
        numerator <- sigma2 * (1 - exp(-gamma2 * l1)) * (1 - exp(-gamma2 * l2)) * exp(-gamma1 * (l1 + l2)/2)
        denominator <- 
          sqrt( ((1 - exp(-gamma2 * l1))^2 + exp(-2 * gamma2 * l1)) * ((1 - exp(-gamma2 * l2))^2 + exp(-2 * gamma2 * l2)) )
        Sigma[l1+1, l2+1] <- numerator/denominator
      }
    }
  }
  return(Sigma)
}

### Function that constructs the data structures for X, Y, and Z
### X is a 3-D array of lags (n_counties x n_days x n_lags). First dim is county, and last 2 define a lag matrix 
### Y is a 2-D matrix (n_counties x n_days) where rows are counties and columns are days
### Z is a 3-D array (n_counties x n_days x # of covariates+1).  First dim is county, last 2 is a covariate matrix for that county (n_days x n_covariates+1)
### D is a 1-D vector of offsets
construct_data_matrices <- function(df, covariates) {
  n_counties <- n_distinct(df$fips)
  n_lags <- sum(grepl('lag_', names(df))) 
  n_days <- 
    df %>% 
    group_by(fips) %>% 
    summarise('n_days' = n()) %>% 
    pull(n_days) %>% 
    max()
  
  n_days <- n_days - (n_lags - 1)
  
  counties <- unique(df$fips)
  lags <- paste0('lag_', 0:(n_lags-1))
  
  ### Initialize storage for X, Y, Z
  X <- 
    array(dim = c(n_counties, n_days, n_lags),
          dimnames = list(counties, 1:n_days, lags))
  
  Y <- matrix(NA, nrow = n_counties, ncol = n_days)
  
  Z <- 
    array(dim = c(n_counties, n_days, length(covariates) + 1),
          dimnames = list(counties, 1:n_days, c('(Intercept)', covariates)))
  
  D <- vector(length = n_counties)
  names(D) <- counties
  
  ### Iterate over counties and add in the X, Y, Z component specific to that county
  for(c in 1:n_counties) {
    X[c,,] <- 
      filter(df, fips == counties[c]) %>% 
      drop_na() %>% 
      select(contains('lag_')) %>% 
      as.matrix()
    

    Y[c,] <- 
      filter(df, fips == counties[c]) %>% 
      drop_na() %>% 
      select(n_cases) %>% 
      as.matrix() %>% 
      t(.)
    
    Z[c,,1] <- 1
    Z[c,,-1] <- 
      filter(df, fips == counties[c]) %>% 
      drop_na() %>% 
      select(all_of(covariates)) %>% 
      as.matrix()
    
    D[c] <- 
      filter(df, fips == counties[c]) %>% 
      drop_na() %>% 
      pull(population) %>% 
      head(1)
    
  }
  
  data <- list('X' = X, 'Y' = Y, 'Z' = Z, 'D' = D)
  return(data)
  
}
