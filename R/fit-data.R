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
library(xtable)
library(ggpubr)
library(tidyselect)

source(here("R", "fit-data-helper-functions.R"))

# for plots
p_cols = c("blue" = "blue4",
            "green" = "forestgreen", 
            "if_bg" = "hotpink1", 
            "if_gb" = "sienna1" , 
            "if_nbg" = "deeppink3", 
            "if_ngb" = "orangered3", 
            "AC" = "deepskyblue", 
            "A-C" = "mediumblue", 
            "-AC" = "darkorange", 
            "-A-C" = "lightcoral")

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
####        Independent trials         ####
#### fit Gaussian P(a,c) - P(a) * P(c) ####
###########################################
# 1. Data
df.ind.trials = data.pe %>% dplyr::select(blue, green, AC, prolific_id, id) %>% 
  mutate(diff = AC - blue * green) %>% 
  filter(str_detect(id, "independent"))

df.ind = bind_rows(df.ind.trials, df.ind.trials %>% mutate(id = "all"))
p.ind = df.ind %>%
  ggplot(aes(color = id)) + geom_density(aes(x=diff)) +
  facet_wrap(~id, ncol = 3, scales = "free") +
  theme(legend.position = "top")


# Fit by using Greta ------------------------------------------------------
# 2. variables and priors
mu = normal(0, 1)
sd = uniform(0, 1)

# likelihood
y <- as_data(df.ind$diff)
distribution(y) <- normal(mu, sd)

# 3. model
m <- model(mu, sd)
draws <- mcmc(m, n_samples = 1000)
mat <- data.frame(matrix(draws[[1]], ncol = 2))
names(mat) <- c("mu", "sd")

# 4. Plots
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
# plot MCMC samples with observed data
x <- seq(-0.25, 0.25, by = 0.01)
y = dnorm(x = x, mean = par.fit$mu_hat, sd = par.fit$sd_hat)
p.ind + geom_line(data = tibble(x=x, y=y), aes(x = x, y = y), color = 'red')
summary(draws)

# plot log-likelihood, posterior predictive
N = df.ind$diff %>% length()
df.ll = mat %>% rowid_to_column("i_sample") %>% group_by(i_sample) %>% 
  mutate(ll.X = dnorm(rnorm(N, mean=mu, sd=sd), mean=mu, sd=sd, log=T) %>% sum())

df.ll.obs = mat %>% rowid_to_column("i_sample") %>% group_by(i_sample) %>% 
  mutate(ll.obs = dnorm(df.ind$diff, mean=mu, sd=sd, log=T) %>% sum())


df.ll %>% ggplot() + geom_density(aes(x=ll.X)) +
  geom_point(data = tibble(ll.obs = mean(df.ll.obs$ll.obs)), 
             color = 'firebrick', size = 2,
             aes(x = ll.obs, y = 0)) +
  labs(x = "log likelihood samples from posterior predictive")


# Fit using webppl --------------------------------------------------------
# fit Gaussian independent-data (P(b,g) - P(b)*P(g))
ind_trials = c("all", df.ind$id %>% unique())
posterior_samples.ind = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  if(trial_id == "all") {
    df.trial = df.ind
  } else {
    df.trial = df.ind %>% filter(id == trial_id)
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
    as_tibble() %>% 
    pivot_wider(names_from = "Parameter", values_from = "value") %>% 
    add_column(id = trial_id)
    # group_by(Parameter, Chain) %>% 
    # mutate(Chain = as.factor(Chain), Parameter = as.character(Parameter))  
  
  return(samples.posterior)
})
save_data(posterior_samples.ind, 
          paste(target_dir, "posterior_samples_ind_diffs.rds", sep=FS))

evs.posterior.ind = posterior_samples.ind %>% 
  pivot_longer(cols = c("mu", "sigma"), names_to = "Parameter", 
               values_to = "value") %>% 
  group_by(id, Parameter) %>%
  summarize(ev = mean(value), .groups = "drop_last") %>% 
  pivot_wider(names_from = "Parameter", values_from = "ev")
  
save_data(evs.posterior.ind, paste(target_dir, "fit_ind_diffs.rds", sep=FS))

# plot data with fitted distributions
x <- seq(-1, 1, by = 0.01)
fitted_vals = evs.posterior.ind %>% 
  mutate(y = list(dnorm(x, mean = mu, sd = sigma)), x = list(x)) %>% 
  unnest(c(x, y))

p.ind + geom_line(data = fitted_vals, aes(x=x, y=y), color = 'firebrick')

# plot priors used: for sd and mu
prior_mu.x = seq(-0.5, 0.5, by = 0.01)
prior_sd.x = seq(params$fit_ind_prior_sigma$a, 
                 params$fit_ind_prior_sigma$b, 
                 by = 0.01)
df.prior_mu = tibble(x = prior_mu.x, 
                     y = dnorm(x, mean = params$fit_ind_prior_mu$mu, 
                               sd = params$fit_ind_prior_mu$sigma), 
                     Parameter = "mu", level = "prior")
df.prior_sd = tibble(x = prior_sd.x, y = 1, Parameter = "sd", level = "prior")

bind_rows(df.prior_mu, df.prior_sd) %>%
  ggplot() + geom_line(aes(x=x, y=y, color = Parameter)) +
  theme(legend.position = "none") + 
  facet_wrap(~Parameter, scales = "free")


# Posterior predictive webppl ---------------------------------------------
samples.pp.ind = group_map(posterior_samples.ind %>% group_by(id), 
                                     function(df.samples, df.grp){
  message(paste(df.grp$id))
  # get empirical data for id-probability combination
  if(df.grp$id == "all") {
    behav.trial = df.ind
  } else {
    behav.trial = df.ind %>% filter(id == df.grp$id)
  }
  
  data_webppl = list(probs = behav.trial$diff, 
                     samples_posterior = list(mu = df.samples$mu,
                                              sigma = df.samples$sigma))
  
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-independent-trials.wppl"),
    data_var = "data",
    data = data_webppl
  )
  
  X_new <- unlist(samples.pp$X_new) %>% matrix(ncol = behav.trial %>% nrow(), byrow=T)
  ll_X_new <- samples.pp$ll_X_new
  ll_X_obs = samples.pp$ll_X_obs
  
  result = tibble(X_new = list(X_new), 
                  ll_X_new = list(ll_X_new), 
                  ll_X_obs = list(ll_X_obs))
  return(result %>% add_column(id = df.grp$id))
}) %>% bind_rows()

# log likelihood plot
ll_plot.ind = plot_pp_ll(samples.pp.ind %>% add_column(p = "")) +
  theme(legend.position = "none")
ll_plot.ind


# Fit independent marginals to Zero-one inflated beta ---------------------
# P(blue), P(green)
probs <- c("blue", "green")
df.ind_marginals = df.ind %>% dplyr::select(blue, green, prolific_id, id)
posterior_samples.ind_marginals = fit_zoib(df.ind_marginals, ind_trials, probs)
evs.ind_marginals = get_evs_zoib_samples(posterior_samples.ind_marginals)

pp.ind_marginals = posterior_predictive_zoib(posterior_samples.ind_marginals, 
                                             df.ind_marginals) 
pp_ll_plot.ind = plot_pp_ll(pp.ind_marginals)
pp_ll_plot.ind + theme(legend.position = "none")

save_data(evs.ind_marginals, paste(target_dir, "fit_ind_marginals.rds", sep=FS))


# Dependent Trials --------------------------------------------------------
################################################
####          Dependent Trials              ####   
#### fit Zero-inflated Beta P(c|a), P(c|¬a) ####
################################################
df.dep = data.pe %>%
  rename(if_bg = `if blue falls green falls`, 
         if_nbg = `if blue does not fall green falls`, 
         if_gb = `if green falls blue falls`, 
         if_ngb = `if green does not fall blue falls`) %>% 
  dplyr::select(prolific_id, id, if_bg, if_nbg, if_gb, if_ngb, blue, green, 
                AC, `A-C`, `-AC`, `-A-C`) %>% 
  filter(!str_detect(id, "independent")) %>% 
  mutate(relation = case_when(str_detect(id, "if2") ~ "if2",
                              str_detect(id, "if1") ~ "if1")) %>% 
  group_by(relation)

# fit Zero-one inflated beta distribution
dep_trials = c("all", "if1", "if2", df.dep$id %>% unique())
probs <- c("if_bg", "if_nbg", "if_gb", "if_ngb", "blue", "green")
posterior_samples.dep = fit_zoib(df.dep, dep_trials, probs)
evs.posterior.dep <- get_evs_zoib_samples(posterior_samples.dep)

# posterior predictive
pp.dep = posterior_predictive_zoib(posterior_samples.dep, df.dep) 
# log likelihood plots
pp_ll_plot.dep0 = plot_pp_ll(pp.dep %>% filter(id %in% c("all", "if1", "if2")))
pp_ll_plot.dep0

pp_ll_plot.dep1 = plot_pp_ll(pp.dep %>% filter(startsWith(id, "if1_")))
pp_ll_plot.dep1

pp_ll_plot.dep2 = plot_pp_ll(pp.dep %>% filter(startsWith(id, "if2_")))
pp_ll_plot.dep2

save_data(evs.posterior.dep, paste(target_dir, "fit_dep.rds", sep=FS))

# for table in thesis
evs.table = evs.posterior.dep %>% 
  dplyr::select(alpha, gamma, shape1, shape2, id, p) %>%
  mutate(across(where(is.numeric), function(x){
   return(round(x, digits = 1))
  }))

evs.table %>% filter(p == "blue")


# Fit single table cells --------------------------------------------------
probs.tbls <- c("AC", "A-C", "-AC", "-A-C")
posterior_samples.table_cells = fit_zoib(df.dep, dep_trials)
save_data(posterior_samples.table_cells, 
          paste(target_dir, "posterior_samples_dep_tbls.rds", sep=FS))

evs.posterior.table_cells <- get_evs_zoib_samples(posterior_samples.table_cells)
save_data(evs.posterior.table_cells, paste(target_dir, "fit_dep_tbls.rds", sep=FS))

# posterior predictive
pp.dep.tbls = posterior_predictive_zoib(posterior_samples.table_cells, df.dep) 
save_data(pp.dep.tbls, paste(target_dir, "pp_samples_dep_tbls.rds", sep=FS))

# log likelihood plots
p_cols.tbls <- p_cols[probs.tbls]
ids <- c("all", "if1", "if2")
pp_ll_plot.dep0 = plot_pp_ll(pp.dep.tbls %>% filter(id %in% ids),
                             p_cols.tbls)
pp_ll_plot.dep0

pp_ll_plot.dep1 = plot_pp_ll(pp.dep.tbls %>% filter(startsWith(id, "if1_")),
                             p_cols.tbls)
pp_ll_plot.dep1

pp_ll_plot.dep2 = plot_pp_ll(pp.dep.tbls %>% filter(startsWith(id, "if2_")), 
                             p_cols.tbls)
pp_ll_plot.dep2

# plot posterior predictive
X_new.if1_uh.ac = pp.dep.tbls %>% 
  filter(id == "if1_uh" & p == "AC") %>% pull(X_new) %>% 
  unlist() %>% matrix(ncol = 89, byrow = T)
X_obs.if1_uh.ac = df.dep %>% filter(id == "if1_uh") %>% pull(AC)
  
ppc_dens_overlay(y=X_obs.if1_uh.ac, yrep = X_new.if1_uh.ac[1:50,])

#########################
####                 ####
#### Some more plots ####
####                 ####
#########################

# Some plots causal power, conditional probabilities, etc. ----------------
# plot rated P(green|blue) vs. P(green|¬blue)
df.dep_cp <-   df.dep %>% 
  mutate(causal_power_bg = (if_bg - if_nbg) / (1 - if_nbg),
         causal_power_gb = (if_gb - if_ngb) / (1 - if_ngb))

df.dep_cp %>% filter(!is.na(if_bg) & !is.na(if_nbg)) %>% 
  ggplot(aes(x = if_bg, y = if_nbg)) + 
  geom_hdr(probs = c(0.99, 0.9, 0.5, 0.25)) +
  facet_wrap(~id) + labs(x = "P(green|blue)", y = "P(green|¬blue)")

df.dep_cp %>% filter(!is.na(if_bg) & !is.na(if_nbg)) %>% 
  mutate(cp_larger_noise = causal_power_bg > if_nbg) %>% 
  group_by(id) %>% summarize(m = mean(cp_larger_noise, na.rm = T))

# look at rated causal powers
df.causal_power = df.dep_cp %>% 
  pivot_longer(cols = c(causal_power_bg, causal_power_gb), 
               names_to = "causal_power", 
               names_prefix = "causal_power_",
               values_to = "theta") %>% 
  mutate(causal_power = case_when(causal_power == "bg" ~ "blue block on green block", 
                                  causal_power == "gb" ~ "green block on blue block")) %>% 
  filter(!is.na(theta) & !is.infinite(theta)) 

df.causal_power %>% 
  group_by(id, causal_power) %>% 
  summarize(val = median(theta), .groups = "drop_last") %>% 
  arrange(desc(val))

df.causal_power %>% filter(causal_power == "blue block on green block") %>% 
  filter(theta>0) %>% 
  ggplot(aes(x=theta, y = if_bg)) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  stat_regline_equation(label.x = 0.3, label.y = 0.15) +
  stat_cor(label.x = 0.3, label.y = 0.05) +
  facet_wrap(~id, ncol=4) +
  labs(x = "causal power blue on green block", y = "P(green | blue)")


df.causal_power %>% filter(causal_power == "blue block on green block") %>% 
  ggplot(aes(x = theta)) + geom_histogram(aes(color = id)) +
  facet_wrap(~id)


