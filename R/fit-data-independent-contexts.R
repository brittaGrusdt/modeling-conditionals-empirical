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

target_dir = paste(here(params$dir_results), "independent-contexts", sep=FS)
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, 
                pe_task.smooth, pe_task, slider) %>% 
  translate_standardized2model() 

# use unsmoothed slider ratings from PE-task!
data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task) %>% 
  pivot_wider(names_from = "utt.standardized", 
              values_from = "pe_task") %>% 
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
df.ind = data.pe %>%
  dplyr::select(blue, green, AC, `A-C`, `-AC`, `-A-C`, prolific_id, id) %>% 
  mutate(diff = AC - blue * green) %>% 
  filter(str_detect(id, "independent"))

#df.ind = bind_rows(df.ind.trials, df.ind.trials %>% mutate(id = "all"))
p.ind = df.ind %>%
  ggplot(aes(color = id)) + geom_density(aes(x=diff)) +
  facet_wrap(~id, ncol = 3, scales = "free") +
  theme(legend.position = "top")

ind_trials = df.ind$id %>% unique()

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
    program_file = here("webppl-model", "posterior-independent-trials-diff.wppl"),
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
save_data(evs.posterior.ind, paste(target_dir, "evs_posterior_ind_diffs.rds", sep=FS))


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


# new 
# posterior 
#trials <- ind_trials[!ind_trials %in% c("independent_ul", "independent_ll")]

# normalize rating if necessary, that is, when AC is (due to rounding) not 
# in necessary possible range, e.g. blue: 0.8, green: 0.8 --> min: 0.6, but
# observed value for AC is 0.59
df.to_normalize = df.ind %>% 
  mutate(min_bg = case_when(blue + green > 1 ~ blue + green - 1,
                            T ~ 0),
         max_bg = pmin(blue, green), 
         diff_ac_min = round(AC, 2) - round(min_bg, 2),
         diff_ac_max = round(AC, 2) - round(max_bg, 2),
         ac_out_range = AC < min_bg | AC > max_bg
         ) %>% 
  # check how often observed is limit of possible interval:
  # filter(ac_out_range & ((diff_ac_min != 0 & diff_ac_max == 0) |
  #                          (diff_ac_min == 0 & diff_ac_max != 0) 
  #                          ))
  # not in range and observed ac is not same as min and max value for ac:
  filter(ac_out_range & diff_ac_min !=0 & diff_ac_max != 0)

df.normalized = df.to_normalize %>% 
  group_by(id, prolific_id) %>% 
  mutate(s = sum(AC, `A-C`, `-AC`, `-A-C`), 
         AC = AC/s,
         `A-C` = `A-C`/s, 
         `-AC` = `-AC`/s, 
         `-A-C` = `-A-C`/s, 
         blue = AC + `A-C`, 
         green = AC + `-AC`, 
         min_bg = case_when(blue + green > 1 ~ blue + green - 1,
                            T ~ 0),
         max_bg = pmin(blue, green),
         ac_out_range = AC < min_bg | AC > max_bg)
# all data with those normalized where necessary
df.behav <- bind_rows(
  anti_join(df.ind, df.normalized %>% dplyr::select(prolific_id, id)),
  df.normalized %>%
    dplyr::select(-s, -min_bg, -max_bg, -ac_out_range, -diff_ac_min, -diff_ac_max)
)

posterior_samples.ind = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  df.trial = df.behav %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-independent.wppl"),
    data_var = "data",
    model_var = "non_normalized_posterior",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
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
  
  return(samples.posterior)
})
pars <- c("alpha_blue", "gamma_blue", "shape1_blue", "shape2_blue", 
          "alpha_green", "gamma_green", "shape1_green", "shape2_green",
          "shape1_bg", "shape2_bg")
          #"ratio_range_variance_bg")
evs.posterior.ind = posterior_samples.ind %>% 
  pivot_longer(cols = all_of(pars), names_to = "Parameter", values_to = "value") %>% 
  group_by(id, Parameter) %>%
  summarize(ev = mean(value), .groups = "drop_last") %>% 
  pivot_wider(names_from = "Parameter", values_from = "ev")

# generate set of independent tables with expected values of posterior distributions
sampled_tables = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  data_webppl <- list(probs = df.behav %>% filter(id == trial_id),
                      evs_params = evs.posterior.ind %>% filter(id == trial_id))
  samples.ind_tables <- webppl(
    program_file = here("webppl-model", "posterior-independent.wppl"),
    data_var = "data",
    model_var = "sample_table",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
    inference_opts = list(method = "forward",
                          samples = 1000)
  ) %>% as_tibble() %>% 
  pivot_wider(names_from = "Parameter", values_from = "value") %>% 
  add_column(id = trial_id)
})

# plot new tables sampled for each independent context with observed data
df.samples <- sampled_tables %>% 
  pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), 
               names_to = "cell", values_to = "value")

df.behav.long = df.behav %>% 
  dplyr::select(AC, `A-C`, `-AC`, `-A-C`, prolific_id, id) %>% 
  pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), 
               names_to = "cell", values_to = "value") 

colors <- c("empirical" = "darkgreen", "sampled" = "firebrick")
plots_new_tables <- map(ind_trials, function(trial_id){
  
  trial.behav <- df.behav.long %>% filter(id == trial_id) %>% 
    dplyr::select(id, cell, value) %>% add_column(data = "empirical")
  trial.samples <- df.samples %>% filter(id == trial_id) %>% 
    dplyr::select(id, cell, value) %>% add_column(data = "sampled")
  df.trial <- bind_rows(trial.behav, trial.samples)
  
  evs.behav <- trial.behav %>% group_by(id, cell) %>% 
    summarize(value = mean(value), .groups = "drop_last") %>% 
    add_column(data = "empirical")
  evs.samples <- trial.samples %>% group_by(id, cell) %>% 
    summarize(value = mean(value), .groups = "drop_last") %>% 
    add_column(data = "sampled")
  df.evs <- bind_rows(evs.behav, evs.samples)
  
  p <-  df.trial %>% 
    ggplot(aes(x=value)) + geom_density(aes(color = data)) + 
    geom_point(data = df.evs, aes(y = 0, color = data), size = 2) +
    facet_wrap(~fct_rev(cell), ncol=2, scales = "free") +
    labs(title = trial_id) +
    scale_color_manual(name = "data", values = colors) + 
    theme(legend.position = "top")
  return(p)
})

plots_new_tables

# posterior predictive 
pp_samples_ll = group_map(posterior_samples.ind %>% group_by(id),
                                  function(df.samples, df.grp){
  message(paste(df.grp$id))
  data_webppl <- list(probs = df.behav %>% filter(id == df.grp$id),
                      samples_posterior = df.samples)
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-independent.wppl"),
    data_var = "data",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS))
  ) 

  result = tibble(ll_X_new = samples.pp$ll_X, 
                  ll_obs = samples.pp$ll_obs, 
                  ll_X_blue = samples.pp$ll_x_blue, 
                  ll_X_green = samples.pp$ll_x_green, 
                  ll_X_bg = samples.pp$ll_x_bg, 
                  ll_obs_blue = samples.pp$ll_obs_blue, 
                  ll_obs_green = samples.pp$ll_obs_green, 
                  ll_obs_bg = samples.pp$ll_obs_bg)
  return(result %>% add_column(id = df.grp$id))
}) %>% bind_rows()

# overall log likelihood of data from posterior predictive
ll_X_obs.mean = pp_samples_ll %>% group_by(id) %>% 
  summarize(ev = mean(ll_obs), .groups = "keep")
p <- pp_samples_ll %>% 
  ggplot(aes(x = ll_X_new)) + geom_density() +
  facet_wrap(~id, scales = "free") +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme(legend.position = "top") +
  labs(x = "log likelihood", y = "density")
p

# log likelihood seperate for P(b), P(g) and P(b,g)
pp_samples_ll.long <- pp_samples_ll %>%
  pivot_longer(cols = c(-id), names_to = "prob", names_prefix = "ll_", 
               values_to = "ll")
# ll pp-samples
df.ll_X = pp_samples_ll.long %>% 
  filter(startsWith(prob, "X_") & prob != "X_new") %>% 
  mutate(prob = str_replace(prob, "X_", ""))

# ev ll observed data given params from posterior
ll_X_obs.mean = pp_samples_ll.long %>% group_by(id, prob) %>% 
  filter(startsWith(prob, "obs_")) %>% 
  summarize(ev = mean(ll), .groups = "drop_last") %>% 
  mutate(prob = str_replace(prob, "obs_", ""))

p <- df.ll_X %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(prob~id, ncol=5, scales = "free") +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme(legend.position = "top") +
  labs(x = "log likelihood", y = "density")
p




# Posterior predictive webppl ---------------------------------------------
fn_posterior_predictive <- function(df.samples, df.grp){
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
    program_file = here("webppl-model", 
                        "posterior-predictive-independent-trials-diff.wppl"),
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
}
pp.ind_diff = group_map(posterior_samples.ind %>% group_by(id), 
                        fn_posterior_predictive) %>% bind_rows()
save_data(pp.ind_diff, paste(target_dir, "pp_samples_ind_diff.rds", sep=FS))

# log likelihood plot
ll_plot.ind = plot_pp_ll(pp.ind_diff %>% add_column(p = ""), c("dummy"="")) +
  theme(legend.position = "none")
ll_plot.ind


# Fit independent marginals to Zero-one inflated beta ---------------------
# P(blue), P(green)
probs <- c("blue", "green", "AC")
df.ind_marginals = df.ind %>% dplyr::select(blue, green, AC, prolific_id, id)
posterior_samples.ind_marginals = fit_zoib(df.ind_marginals, ind_trials, probs)
evs.ind_marginals = get_evs_zoib_samples(posterior_samples.ind_marginals)
pp.ind_marginals = posterior_predictive_zoib(posterior_samples.ind_marginals, 
                                             df.ind_marginals) 
save_data(posterior_samples.ind_marginals, 
          paste(target_dir, "posterior_samples_ind_marginals.rds", sep=FS))
save_data(evs.ind_marginals, 
          paste(target_dir, "evs_posterior_ind_marginals.rds", sep=FS))
save_data(pp.ind_marginals, paste(target_dir, "pp_samples_ind_marginals.rds", sep=FS))

pp_ll_plot.ind = plot_pp_ll(pp.ind_marginals %>% filter(id!="all"), p_cols[probs])
pp_ll_plot.ind + theme(legend.position = "none")


# Fit independent tables by cells zero-one inflated beta ------------------
probs.tbls <- c("AC", "A-C", "-AC", "-A-C")
df.ind_tbls = data.pe %>% dplyr::select(AC, `A-C`, `-AC`, `-A-C`, prolific_id, id)
posterior_samples.ind_tbls = fit_zoib(df.ind_tbls, ind_trials, probs.tbls)
evs.ind_tbls = get_evs_zoib_samples(posterior_samples.ind_tbls)
pp.ind_tbls = posterior_predictive_zoib(posterior_samples.ind_tbls, df.ind_tbls) 

save_data(posterior_samples.ind_tbls, 
          paste(target_dir, "posterior_samples_ind_tbls.rds", sep=FS))
save_data(evs.ind_tbls, paste(target_dir, "evs_posterior_ind_tbls.rds", sep=FS))
save_data(pp.ind_tbls, paste(target_dir, "pp_samples_ind_tbls.rds", sep=FS))


pp_ll_plot.ind_tbls = plot_pp_ll(
  pp.ind_tbls %>% 
    mutate(p = factor(p, levels = c('AC', 'A-C', '-AC', '-A-C'))) %>% 
    filter(id != "all"), 
  p_cols[probs.tbls]
)
pp_ll_plot.ind_tbls

context <- "independent_ll"
prob <- "-A-C"
y_rep = pp.ind_tbls %>% filter(id == !!context & p == !!prob) %>% pull(X_new) %>% 
  unlist() %>% matrix(ncol = data.pe %>% filter(id == !!context) %>%
                        nrow(), byrow=T)
ppc_dens_overlay(data.pe %>% filter(id == !!context) %>% pull(prob), 
                 y_rep[1:50, ])

# # summed log-likelihoods across table cells
# pp.across_cells.obs = pp.ind_tbls %>% 
#   dplyr::select(id, ll_X_obs, p) %>% 
#   pivot_wider(names_from = "p", values_from = "ll_X_obs") %>%
#   unnest(c(`AC`, `A-C`, `-AC`, `-A-C`)) %>% 
#   mutate(ll_X_obs = `AC` + `A-C` + `-AC` + `-A-C`) %>% 
#   dplyr::select(id, ll_X_obs) %>% arrange(id)
# 
# pp.across_cells.new = pp.ind_tbls %>% 
#   dplyr::select(id, ll_X_new, p) %>% 
#   pivot_wider(names_from = "p", values_from = "ll_X_new") %>%
#   unnest(c(`AC`, `A-C`, `-AC`, `-A-C`)) %>% 
#   mutate(ll_X_new = `AC` + `A-C` + `-AC` + `-A-C`) %>% 
#   dplyr::select(id, ll_X_new) %>% arrange(id)
# 
# pp.across_cells = bind_cols(pp.across_cells.obs %>% dplyr::select(-id), 
#                             pp.across_cells.new) %>% 
#   group_by(id) %>% summarize(ll_X_obs = list(ll_X_obs), ll_X_new = list(ll_X_new)) 
# 
# pp.ind_across_cells = left_join(
#   pp.across_cells, 
#   pp.ind_tbls %>% dplyr::select(X_new, id) %>% distinct()
# )
#   
# pp_ll_plot.ind_across_cells = plot_pp_ll(
#   pp.ind_across_cells %>% mutate(p = "") %>% filter(id != "all"),
#   c("dummy"="")
# )
# pp_ll_plot.ind_across_cells + theme(legend.position = "none")
# 


# Which model is a better fit for the data? -------------------------------

#' ll for zoib data params = evs posterior
get_log_likelihood_data = function(df, evs, probs){

  df <- left_join(
    df %>% 
      pivot_longer(cols = probs, names_to = "p", values_to = "p_hat"), 
    evs %>% filter(id != "all"), 
    by = c("id", "p")
  ) 
  
  df.zero_one = df %>% filter(p_hat %in% c(1,0)) %>%
    mutate(ll = dbinom(p_hat, 1, alpha, log=T) + dbinom(p_hat, 1, gamma, log=T)) %>% 
    dplyr::select(id, prolific_id, p, ll)
  
  df.interval = df %>% filter(!p_hat %in% c(1,0)) %>% 
    mutate(ll = dbeta(p_hat, shape1 = shape1, shape2 = shape2, log=T)) %>% 
    dplyr::select(id, prolific_id, p, ll)
  
  df.ll = bind_rows(df.zero_one, df.interval) %>% 
    group_by(id, prolific_id) %>% 
    summarize(ll = sum(ll), .groups = "drop_last") %>% filter(id != "all")
  
  df.ll.ids = df.ll %>% summarize(ll = sum(ll))
  return(df.ll.ids) 
}

# log likelihood marginals, P(a,c)-P(a)*P(c)
ll_marginals = get_log_likelihood_data(df.ind, evs.ind_marginals %>% filter(id!="all"), 
                                       c("blue", "green")) %>% 
  rename(ll.marginals = ll)

ll_marginals_diff = left_join(
  ll_marginals, 
  left_join(df.ind %>% dplyr::select(prolific_id, id, diff), 
  evs.posterior.ind) %>% 
    mutate(ll = dnorm(diff, mean=mu, sd = sigma, log = T)) %>% 
    group_by(id) %>% summarize(ll.diff = sum(ll))
) %>% mutate(ll.marginal_diff = ll.marginals + ll.diff)


df.tbls <- data.pe %>% dplyr::select(AC, `A-C`, `-AC`, `-A-C`, id, prolific_id) %>% 
  filter(startsWith(id, "independent"))
ll_tbls = get_log_likelihood_data(df.tbls, evs.ind_tbls %>% filter(id!="all"),
                                  c("AC", "A-C", "-AC", "-A-C")) %>% 
  rename(ll.tbls = ll)


bf = left_join(ll_tbls, ll_marginals_diff, by = "id") %>% 
  mutate(bf = exp(ll.tbls)/(exp(ll.marginal_diff)))

