library(tidyverse)
library(furrr)

### Set up Parallelization
n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 10 * 1024^3)

### Ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 14),
                  plot.caption = element_text(size = 10),
                  legend.position = "bottom"))

cvd_draws <- read_rds('~/Dropbox (Harvard University)/Bayesian/draws_cvd.rds')
resp_draws  <- read_rds('~/Dropbox (Harvard University)/Bayesian/draws_resp.rds')


### CVD Analysis
### (1) National Averages
cvd_mu <- 
  cvd_draws$mu %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'lags',
               values_to = 'mu',
               names_prefix = 'lag_') %>% 
  group_by('lag' = as.numeric(lags)) %>% 
  summarise('mean' = mean(exp(mu * 10) - 1),
            'lower' = quantile(exp(mu * 10) - 1, 0.025),
            'upper' = quantile(exp(mu * 10) - 1, 0.975)) 

ggplot(cvd_mu, aes(x = lag, y = mean)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = 0:7) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = '# of Lags',
       y = expression(paste('% Increase in Admissions due to 10', mu, 'g/m'^3, ' increase in ', PM[2.5])),
       title = 'National Average Distributed Lag Function',
       subtitle = 'Coronary Artery Disease') + 
  theme(legend.position = 'none')

ggsave('figures/cad_national_avg.png', height = 9/1.2, width = 16/1.2)

cvd_draws$mu %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'lags',
               values_to = 'mu') %>% 
  ggplot(aes(x = iteration, y = mu)) + 
  facet_wrap(~lags, ncol = 4) + 
  geom_line(aes(col = lags)) + 
  labs(x = 'Iteration',
       y = 'mu',
       title = 'Trace Plots for National Avg. Distributed Lag',
       subtitle = 'Coronary Artery Disease')

ggsave('figures/cad_trace_plots.png', height = 9/1.2, width = 16/1.2)

       

### (2) County Averages
cvd_theta <- 
  future_map_dfr(dimnames(cvd_draws$theta)[[2]], ~{
    cvd_draws$theta[,.x,] %>% 
      as_tibble() %>% 
      mutate('iteration' = 1:nrow(.)) %>% 
      pivot_longer(cols = -iteration,
                   names_to = 'lags',
                   values_to = 'mu',
                   names_prefix = 'lag_') %>% 
      group_by('lag' = as.numeric(lags)) %>% 
      summarise('mean' = mean(exp(mu * 10) - 1),
                'lower' = quantile(exp(mu * 10) - 1, 0.025),
                'upper' = quantile(exp(mu * 10) - 1, 0.975)) %>% 
      mutate('fips' = .x)
  })



ggplot(cvd_theta, aes(x = lag, y = mean)) +
  facet_wrap(~fips) + 
  geom_line() +
  geom_hline(yintercept = 0, lty = 2, col = 'red') + 
  geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = 0:7) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = '# of Lags',
       y = expression(paste('% Increase in Admissions due to 10', mu, 'g/m'^3, ' increase in ', PM[2.5])),
       title = 'County Average Distributed Lag Function',
       subtitle = 'Coronary Artery Disease')

ggsave('figures/cad_county_avg.png', height = 9 * 1.2, width = 16 * 1.2)

### (3) Gamma and Eta
cvd_draws$gamma %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'term',
               values_to = 'estimate') %>% 
  bind_rows(
    cvd_draws$eta %>% 
      as_tibble() %>% 
      mutate('iteration' = 1:nrow(.)) %>% 
      pivot_longer(cols = -iteration,
                   names_to = 'term',
                   values_to = 'estimate') 
  ) %>% 
  ggplot(aes(x = iteration, y = estimate)) + 
  facet_wrap(~term, scales = 'free') + 
  geom_line(aes(col = term)) + 
  labs(x = 'Iteration',
       y = 'Estimate',
       title = 'Trace Plots for Gamma and Eta',
       subtitle = 'Coronary Artery Disease') + 
  theme(legend.position = 'none')

ggsave('figures/cad_eta_gamma.png', height = 9/1.2, width = 16/1.2)


### RESP Analysis
### (1) National Averages
resp_mu <- 
  resp_draws$mu %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'lags',
               values_to = 'mu',
               names_prefix = 'lag_') %>% 
  group_by('lag' = as.numeric(lags)) %>% 
  summarise('mean' = mean(exp(mu * 10) - 1),
            'lower' = quantile(exp(mu * 10) - 1, 0.025),
            'upper' = quantile(exp(mu * 10) - 1, 0.975)) 

ggplot(resp_mu, aes(x = lag, y = mean)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = 0:7) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = '# of Lags',
       y = expression(paste('% Increase in Admissions due to 10', mu, 'g/m'^3, ' increase in ', PM[2.5])),
       title = 'National Average Distributed Lag Function',
       subtitle = 'Chronic Obstructive Pulmonary Disease')

ggsave('figures/copd_national_avg.png', height = 9/1.2, width = 16/1.2)

resp_draws$mu %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'lags',
               values_to = 'mu') %>% 
  ggplot(aes(x = iteration, y = mu)) + 
  facet_wrap(~lags, ncol = 4) + 
  geom_line(aes(col = lags)) + 
  labs(x = 'Iteration',
       y = 'mu',
       title = 'Trace Plots for National Avg. Distributed Lag',
       subtitle = 'Chronic Obstructive Pulmonary Disease')

ggsave('figures/resp_trace_plots.png', height = 9/1.2, width = 16/1.2)

### (2) County Averages
resp_theta <- 
  future_map_dfr(dimnames(resp_draws$theta)[[2]], ~{
    cvd_draws$theta[,.x,] %>% 
      as_tibble() %>% 
      mutate('iteration' = 1:nrow(.)) %>% 
      pivot_longer(cols = -iteration,
                   names_to = 'lags',
                   values_to = 'mu',
                   names_prefix = 'lag_') %>% 
      group_by('lag' = as.numeric(lags)) %>% 
      summarise('mean' = mean(exp(mu * 10) - 1),
                'lower' = quantile(exp(mu * 10) - 1, 0.025),
                'upper' = quantile(exp(mu * 10) - 1, 0.975)) %>% 
      mutate('fips' = .x)
  })

ggplot(resp_theta, aes(x = lag, y = mean)) +
  facet_wrap(~fips) + 
  geom_line() +
  geom_hline(yintercept = 0, lty = 2, col = 'red') + 
  geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = 0:7) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = '# of Lags',
       y = expression(paste('% Increase in Admissions due to 10', mu, 'g/m'^3, ' increase in ', PM[2.5])),
       title = 'County Average Distributed Lag Function',
       subtitle = 'Chronic Obstructive Pulmonary Disease')
ggsave('figures/copd_county_avg.png', height = 9 * 1.2, width = 16 * 1.2)

### (3) Gamma and eta
resp_draws$gamma %>% 
  as_tibble() %>% 
  mutate('iteration' = 1:nrow(.)) %>% 
  pivot_longer(cols = -iteration,
               names_to = 'term',
               values_to = 'estimate') %>% 
  bind_rows(
    resp_draws$eta %>% 
      as_tibble() %>% 
      mutate('iteration' = 1:nrow(.)) %>% 
      pivot_longer(cols = -iteration,
                   names_to = 'term',
                   values_to = 'estimate') 
  ) %>% 
  ggplot(aes(x = iteration, y = estimate)) + 
  facet_wrap(~term, scales = 'free') + 
  geom_line(aes(col = term)) + 
  labs(x = 'Iteration',
       y = 'Estimate',
       title = 'Trace Plots for Gamma and Eta',
       subtitle = 'Chronic Obstructive Pulmonary Disease') + 
  theme(legend.position = 'none')

ggsave('figures/copd_eta_gamma.png', height = 9/1.2, width = 16/1.2)



