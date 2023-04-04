library(tidyverse)


source('helpers.R')
source('sampler.R')

### Read in Data 
### First 2 columns of cvd.csv and resp.csv are just row #s so get rid of those

set.seed(6969)
cases <- read_rds('sample_data/cases.rds')
cases[[1]] <- sample(cases[[1]])
cases[[2]] <- sample(cases[[2]])

### Some coeff to generate fake data
### This will be useful for debugging so we know what the truth is 
theta_sim <- 0.0002 * seq(1, 0, length.out = 8)
beta_sim <- c(-10, 0.005, 0.002)


cvd <- 
  read_csv('sample_data/cvd.csv') %>% 
  select(-c(1:2)) 

resp <- 
  read_csv('sample_data/resp.csv') %>% 
  select(-c(1:2)) 

load('sample_data/gridmet_final.Rdata')

confounders <- 
  inner_join(bind_rows(gridMET_subset[[1]], gridMET_subset[[2]]),
             bind_rows(gridMET_subset[[3]], gridMET_subset[[4]]),
             by = c('date', 'fips')) %>% 
  mutate('date' = as.Date(date))

### Impute confounders by mean of measurements on that day
confounders <-
  confounders %>% 
  group_by(date) %>% 
  mutate('mean_rmax' = mean(rmax, na.rm = T),
         'mean_tmmx' = mean(tmmx, na.rm = T)) %>% 
  ungroup() %>% 
  mutate('rmax' = ifelse(is.na(rmax), mean_rmax, rmax),
         'tmmx' = ifelse(is.na(tmmx), mean_tmmx, tmmx))

### Split our df into one fips at a time and then map our function to create the 
### lag df over each fips, and bind them back together via bind_rows
### (this is what map_dfr does) 
cvd <- 
  cvd %>% 
  group_by(fips) %>% 
  group_split() %>% 
  map_dfr(~create_lag_df(.x, n_lags = 7, lag_var = 'PM25_mean_pred')) %>% 
  mutate('date' = as.Date(paste(year, month, day), '%Y %m %d')) %>% 
  inner_join(confounders, by = c('date', 'fips')) 

### Add in some Fake Cases
M <- 
  cvd %>% 
  mutate('intercept' = 1) %>% 
  select(contains('lag'), 'intercept', 'tmmx', 'rmax') %>% 
  as.matrix() 
lambda <- exp(M %*% c(theta_sim, beta_sim) + log(cvd$population))
cvd$n_cases <- rpois(n = nrow(cvd), lambda = lambda)

resp <- 
  resp %>% 
  group_by(fips) %>% 
  group_split() %>% 
  map_dfr(~create_lag_df(.x, n_lags = 7, lag_var = 'PM25_mean_pred')) %>% 
  mutate('date' = as.Date(paste(year, month, day), '%Y %m %d')) %>% 
  inner_join(confounders, by = c('date', 'fips'))

### Add in some Fake Cases
M <- 
  resp %>% 
  mutate('intercept' = 1) %>% 
  select(contains('lag'), 'intercept', 'tmmx', 'rmax') %>% 
  as.matrix() 
lambda <- exp(M %*% c(theta_sim, beta_sim) + log(cvd$population))
resp$n_cases <- rpois(n = nrow(resp), lambda = lambda)

### Create X, Y, and Z matrices for each county
cvd_data <- construct_data_matrices(df = cvd, covariates = c('tmmx', 'rmax'))
resp_data <- construct_data_matrices(df = resp, covariates = c('tmmx', 'rmax'))



X <- cvd_data$X
Y <- cvd_data$Y
Z <- cvd_data$Z
D <- cvd_data$D

t0 <- Sys.time()
draws <- 
  gibbs_sampler(X = cvd_data$X, 
                Y = cvd_data$Y,
                Z = cvd_data$Z, 
                D = cvd_data$D,
                burn_in = 10, 
                n_samples = 500)
t1 <- Sys.time()
t1 - t0
