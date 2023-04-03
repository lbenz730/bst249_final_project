library(tidyverse)


source('helpers.R')
source('sampler.R')

### Read in Data 
### First 2 columns of cvd.csv and resp.csv are just row #s so get rid of those
### Cases[[1]] is fake outcome for CVD. I'll randomly shuffle it 
### Cases[[2]] is fake outcome for resp. I'll randomly shuffle it 
set.seed(6969)
cases <- read_rds('sample_data/cases.rds')
cases[[1]] <- sample(cases[[1]])
cases[[2]] <- sample(cases[[2]])

cvd <- 
  read_csv('sample_data/cvd.csv') %>% 
  select(-c(1:2)) %>% 
  mutate('n_cases' = cases[[1]])

resp <- 
  read_csv('sample_data/resp.csv') %>% 
  select(-c(1:2)) %>% 
  mutate('n_cases' = cases[[2]])


### Split our df into one fips at a time and then map our function to create the 
### lag df over each fips, and bind them back together via bind_rows
### (this is what map_dfr does) 
cvd <- 
  cvd %>% 
  group_by(fips) %>% 
  group_split() %>% 
  map_dfr(~create_lag_df(.x, n_lags = 7, lag_var = 'PM25_mean_pred'))

resp <- 
  resp %>% 
  group_by(fips) %>% 
  group_split() %>% 
  map_dfr(~create_lag_df(.x, n_lags = 7, lag_var = 'PM25_mean_pred'))

### For now make some fake noise for Z
cvd$z1 <- runif(nrow(cvd), max = 100)
cvd$z2 <- runif(nrow(cvd), max = 100)
resp$z1 <- runif(nrow(resp), max = 100)
resp$z2 <- runif(nrow(resp), max = 100)


### Create X, Y, and Z matrices for each county
cvd_data <- construct_data_matrices(df = cvd, covariates = c('z1', 'z2'))
resp_data <- construct_data_matrices(df = resp, covariates = c('z1', 'z2'))



X <- cvd_data$X
Y <- cvd_data$Y
Z <- cvd_data$Z
D <- cvd_data$D

