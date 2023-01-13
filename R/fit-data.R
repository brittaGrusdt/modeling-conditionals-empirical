library(brms)
library(greta)
library(here)
library(tidyverse)
library(tibble)
library(broom)
library(broom.mixed)
library(knitr)
library(ExpDataWrangling)
library(ModelUtils)
library(rwebppl)
library(bayesplot)
# Data --------------------------------------------------------------------
active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = here(params$dir_results)
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
  translate_standardized2model() 

data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task) %>% 
  pivot_wider(names_from = "utt.standardized", values_from = "pe_task") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`) %>% 
  get_controlled_factors() %>% dplyr::select(-relation_type) %>% 
  rename(blue = `blue falls`, green = `green falls`) %>% 
  mutate(blue_zero_one = (blue == 0 | blue == 1),
         green_zero_one = (green == 0 | green == 1))



# Independent Trials ------------------------------------------------------
###########################################
#### fit Gaussian P(a,c) - P(a) * P(c) ####
###########################################
# 1. Data
df.gaussian = data.pe %>% dplyr::select(blue, green, AC, prolific_id, id) %>% 
  mutate(diff = AC - blue * green) %>% 
  filter(str_detect(id, "independent"))
p.ind = df.gaussian %>% ggplot() + geom_density(aes(x=diff))

# 2. variables and priors
mu = normal(0, 1)
sd = uniform(0, 1)

# likelihood
y <- as_data(df.gaussian$diff)
distribution(y) <- normal(mu, sd)

# 3. model
m <- model(mu, sd)
draws <- mcmc(m, n_samples = 1000)
mat <- data.frame(matrix(draws[[1]], ncol = 2))
names(mat) <- c("mu", "sd")

# plot posterior P(mu|X)
p.mu <- ggplot(mat, aes(x=mu)) + 
  geom_histogram(aes(y=..density..), binwidth=.0005, colour="black", 
                 fill="white") +
  geom_density(alpha = .2, fill = '#FF6666')
p.mu
# plot posterior P(sd|X)
p.sd <- ggplot(mat, aes(x=sd)) + 
  geom_histogram(aes(y=..density..), binwidth=.0005, colour="black", 
                 fill="white") +
  geom_density(alpha = .2, fill = '#FF6666')
p.sd
fit <- opt(m)
par.fit = tibble(mu_hat = fit$par$mu, 
                 sd_hat = fit$par$sd)
par.fit
x <- seq(-0.25, 0.25, by = 0.01)
y = dnorm(x = x, mean = par.fit$mu_hat, sd = par.fit$sd_hat)
p.ind + geom_line(data = tibble(x=x, y=y), aes(x = x, y = y), color = 'red')
summary(draws)

# log-likelihood plots posterior predictive
N = df.gaussian$diff %>% length()
df.ll = mat %>% rowid_to_column("i_sample") %>% group_by(i_sample) %>% 
  mutate(ll.X = dnorm(rnorm(N, mean=mu, sd=sd), mean=mu, sd=sd, log=T) %>% sum())

df.ll.obs = mat %>% rowid_to_column("i_sample") %>% group_by(i_sample) %>% 
  mutate(ll.obs = dnorm(df.gaussian$diff, mean=mu, sd=sd, log=T) %>% sum())


df.ll %>% ggplot() + geom_density(aes(x=ll.X)) +
  geom_point(data = tibble(ll.obs = mean(df.ll.obs$ll.obs)), 
             color = 'firebrick', size = 2,
             aes(x = ll.obs, y = 0)) +
  labs(x = "log likelihood samples from posterior predictive")


# Same with webppl --------------------------------------------------------
# fit Gaussian independent-data
ind_trials = c("all", df.gaussian$id %>% unique())
evs.posterior = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  if(trial_id == "all") {
    df.trial = df.gaussian
  } else {
    df.trial = df.gaussian %>% filter(id == trial_id)
  }
  data_webppl = list(probs = df.trial$diff, 
                     prior_mu = params$fit_ind_prior_mu, 
                     prior_sigma = params$fit_ind_prior_sigma)
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-independent-trials.wppl"),
    data_var = "data",
    model_var = "non_normalized_posterior",
    data = data_webppl,
    inference_opts = list(method = "MCMC",
                          samples = 1000,
                          lag = 2,
                          burn = 10000,
                          verbose = T),
    chains = 4, 
    cores = 4) %>%
    as_tibble() %>% group_by(Parameter, Chain) %>% 
    mutate(Chain = as.factor(Chain), Parameter = as.character(Parameter))  
  
  evs.posterior = samples.posterior %>% group_by(Parameter) %>% 
    summarize(ev = mean(value)) %>% 
    pivot_wider(names_from = "Parameter", values_from = "ev") %>% 
    add_column(id = trial_id)
  return(evs.posterior)
})
evs.posterior

x <- seq(-1, 1, by = 0.01)
y = dnorm(x = x, mean = evs.posterior$mu, sd = evs.posterior$sigma)
p.ind + geom_line(data = tibble(x=x, y=y), aes(x = x, y = y), color = 'red')


# plot priors used: for sd and mu
x.prior_mu = seq(-0.5, 0.5, by = 0.01)
x.prior_sd = seq(data_webppl$prior_sigma$a, 
                 data_webppl$prior_sigma$b, 
                 by = 0.01)
df.prior_mu = tibble(x = x.prior_mu, 
                     y = dnorm(x, mean = data_webppl$prior_mu$mu, 
                               sd = data_webppl$prior_mu$sigma), 
                     Parameter = "mu", level = "prior")
df.prior_sd = tibble(x = x.prior_sd, y = 1, Parameter = "sd", level = "prior")

bind_rows(df.prior_mu, df.prior_sd) %>% 
  ggplot() + geom_line(aes(x=x, y=y, color = Parameter)) +
  geom_density(data = posterior_samples, aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") 


# Posterior predictive webppl ---------------------------------------------
samples.mu = samples.posterior %>% filter(Parameter == "mu") %>% pull(value)
samples.sigma = samples.posterior %>% filter(Parameter == "sigma") %>% pull(value)

data_webppl = list(probs = df.gaussian$diff, 
                   samples_posterior = list(mu = samples.mu,
                                            sigma = samples.sigma))
samples.pp <- webppl(
  program_file = here("webppl-model", "posterior-predictive-independent-trials.wppl"),
  data_var = "data",
  data = data_webppl
) %>% as.matrix()  

ppc_dens_overlay(data_webppl$probs, samples.pp[1:50,])

# log likelihood of samples from posterior predictive (return ll_X_new by wppl)
tibble(ll_X = samples.pp) %>% ggplot(aes(x = ll_X)) + geom_density()

# Dependent Trials --------------------------------------------------------

###############################################
#### fit Zero-inflated Beta P(c|a), P(c|¬a) ###
###############################################
df.dep = data.pe %>%
  rename(if_bg = `if blue falls green falls`, 
         if_nbg = `if blue does not fall green falls`, 
         if_gb = `if green falls blue falls`, 
         if_ngb = `if green does not fall blue falls`) %>% 
  dplyr::select(prolific_id, id, if_bg, if_nbg, if_gb, if_ngb, blue, green) %>% 
  filter(!str_detect(id, "independent")) %>% 
  mutate(relation = case_when(str_detect(id, "if2") ~ "if2",
                              str_detect(id, "if1") ~ "if1")) %>% 
  group_by(relation)

# plot rated P(green|blue) vs. P(green|¬blue)
df.dep %>% filter(!is.na(if_bg) & !is.na(if_nbg)) %>% 
  ggplot(aes(x = if_bg, y = if_nbg)) + 
  geom_hdr(probs = c(0.99, 0.9, 0.5, 0.25)) +
  facet_wrap(~id) + labs(x = "P(green|blue)", y = "P(green|¬blue)")
  

dep_trials = c("all", "if1", "if2", df.dep$id %>% unique())
evs.posterior.dep = map_dfr(dep_trials, function(trial_id){
  message(trial_id)
  if(trial_id == "all") {
    df.trial = df.dep
  } else if(trial_id == "if1") {
    df.trial = df.dep %>% filter(relation == "if1")
  } else if(trial_id == "if2") {
    df.trial = df.dep %>% filter(relation == "if2")
  } else {
    df.trial = df.dep %>% filter(id == trial_id)
  }
  probs <- c("if_bg", "if_nbg", "if_gb", "if_ngb")
  evs.posterior_trial = map_dfr(probs, function(p) {
    df.trial_p <- df.trial %>% rename(p = !!p) %>% filter(!is.na(p))
    data_webppl = list(probs = df.trial_p)
    samples.posterior <- webppl(
      program_file = here("webppl-model", "posterior-dependent-trials.wppl"),
      data_var = "data",
      model_var = "non_normalized_posterior",
      data = data_webppl,
      inference_opts = list(method = "MCMC",
                            samples = 1000,
                            lag = 2,
                            burn = 10000,
                            verbose = T),
      chains = 4, 
      cores = 4) %>% as_tibble() 
    
    evs.posterior = samples.posterior %>% group_by(Parameter) %>% 
      summarize(ev = mean(value)) %>% 
      pivot_wider(names_from = "Parameter", values_from = "ev") %>% 
      add_column(id = trial_id, p = p)
    return(evs.posterior)
  })
  return(evs.posterior_trial)
})
evs.posterior.dep
save_data(evs.posterior.dep, paste(target_dir, "fit_dep.rds", sep=""))


# Posterior predictive dependent trials -----------------------------------
# (for one single run, prob: P(g|b), if1_uh)
trial_id = "if_uh" #(run once inner loop above)
p <- "if_bg"
samples.shape1 = samples.posterior %>% filter(Parameter == "shape1") %>% pull(value)
samples.shape2 = samples.posterior %>% filter(Parameter == "shape2") %>% pull(value)
samples.alpha = samples.posterior %>% filter(Parameter == "alpha") %>% pull(value)
samples.gamma = samples.posterior %>% filter(Parameter == "gamma") %>% pull(value)

data_webppl = list(probs = df.trial[[p]], 
                   samples_posterior = list(shape1 = samples.shape1,
                                            shape2 = samples.shape2,
                                            alpha = samples.alpha,
                                            gamma = samples.gamma))
samples.pp <- webppl(
  program_file = here("webppl-model", "posterior-predictive-dependent-trials.wppl"),
  data_var = "data",
  data = data_webppl
) %>% as.matrix()

ppc_dens_overlay(df.trial[[p]], samples.pp[1:50,])

# log likelihood of samples from posterior predictive (return ll_X_new by wppl)
tibble(ll_X = samples.pp) %>% ggplot(aes(x = ll_X)) + geom_density()

#expected log likelihood for draws from posterior for observed data 
# (other return value webppl)
mean(samples.pp) 




