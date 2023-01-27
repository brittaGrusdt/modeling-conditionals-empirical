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
# (data model with zero-one inflated beta!)
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

df.ind = data.pe %>%
  dplyr::select(blue, green, AC, `A-C`, `-AC`, `-A-C`, prolific_id, id) %>% 
  mutate(diff = AC - blue * green) %>% 
  filter(str_detect(id, "independent"))

ind_trials = df.ind$id %>% unique()

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




# Posterior ---------------------------------------------------------------
posterior_samples.ind = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  df.trial = df.behav %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-independent.wppl"),
    data_var = "data",
    model_var = "non_normalized_posterior",
    random_seed = params$seed_webppl,
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

save_data(evs.posterior.ind, 
          paste(target_dir, "evs-posterior-independent-data.rds", sep=FS))

# Generate new independent tables -----------------------------------------
# generate set of independent tables with expected values of posterior distributions
sampled_tables = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  data_webppl <- list(probs = df.behav %>% filter(id == trial_id),
                      evs_params = evs.posterior.ind %>% filter(id == trial_id))
  samples.ind_tables <- webppl(
    program_file = here("webppl-model", "posterior-independent.wppl"),
    random_seed = params$seed_webppl,
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
    theme_minimal() +
    theme(legend.position = "top")
  return(p)
})

plots_new_tables


# Posterior predictive ----------------------------------------------------
# Log likelihood plots
pp_samples_ll = group_map(posterior_samples.ind %>% group_by(id),
                          function(df.samples, df.grp){
                            message(paste(df.grp$id))
                            data_webppl <- list(probs = df.behav %>% filter(id == df.grp$id),
                                                samples_posterior = df.samples)
                            samples.pp <- webppl(
                              program_file = here("webppl-model", "posterior-predictive-independent.wppl"),
                              random_seed = params$seed_webppl,
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
  theme_minimal() +
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
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "log likelihood", y = "density")
p
