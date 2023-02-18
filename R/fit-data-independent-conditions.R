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
library(ggthemes)

source(here("R", "helpers-data-models.R"))

# for plots
theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
prob_names <- c("blue"="P(b)", "green" = "P(g)")
trial_names <- c("independent_hh" = "ind:HH", 
                 "independent_ll" = "ind:LL", 
                 "independent_uh" = "ind:UH", 
                 "independent_ul" = "ind:UL", 
                 "independent_hl" = "ind:HL")

# Data --------------------------------------------------------------------
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "independent-contexts", 
                   "zoib-model", sep=FS)
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
  mutate(bg_ind = blue * green,
         diff = AC - bg_ind,
         min_bg = case_when(blue + green > 1 ~ blue + green - 1, T ~ 0),
         max_bg = pmin(blue, green), 
         min_is_max_bg = round(min_bg, 2) == round(max_bg, 2),
         diff_ac_min = round(AC, 2) - round(min_bg, 2),
         diff_ac_max = round(max_bg, 2) - round(AC, 2)
         ) %>% 
  filter(str_detect(id, "independent"))

ind_trials = df.ind$id %>% unique()

# normalize rating if necessary, that is, when AC is (due to rounding) not 
# in necessary possible range, e.g. blue: 0.8, green: 0.8 --> min: 0.6, but
# observed value for AC is 0.59
df.to_normalize = df.ind %>% 
  mutate(ac_out_range = AC < min_bg | AC > max_bg) %>% 
  filter(ac_out_range & !min_is_max_bg & diff_ac_min !=0 & diff_ac_max != 0)
  # check how often observed is limit of possible interval:
  # filter(diff_ac_max == 0 | diff_ac_min == 0)
  # not in range and observed ac is not same as min and max value for ac:

df.normalized = df.to_normalize %>% 
  group_by(id, prolific_id) %>% 
  mutate(s = sum(AC, `A-C`, `-AC`, `-A-C`), 
         AC = AC/s,
         `A-C` = `A-C`/s, 
         `-AC` = `-AC`/s, 
         `-A-C` = `-A-C`/s, 
         blue = AC + `A-C`, 
         green = AC + `-AC`, 
         bg_ind = blue * green,
         diff = AC - bg_ind,
         min_bg = case_when(blue + green > 1 ~ blue + green - 1, T ~ 0),
         max_bg = pmin(blue, green)
  )

# all data with those normalized where necessary
df.behav <- bind_rows(
  anti_join(df.ind, df.normalized %>% dplyr::select(prolific_id, id)),
  df.normalized %>%
    dplyr::select(-s, -ac_out_range)
)

# check ratio difference P(b,g)-P(b)*P(g)
df.behav %>% 
  mutate(
    max_sub = bg_ind - min_bg, # at most this subtracted from perfect indep P(a,c)=P(b)*P(g)
    max_add = max_bg - bg_ind, # at most this added to perfect indep P(a,c)=P(b)*P(g)
    delta = case_when(diff < 0 ~ -diff/max_sub,
                      diff > 0 ~ diff/max_add,
                      T ~ 0)
    ) %>% 
  ggplot(aes(x = delta)) + geom_density() +
  facet_wrap(~id, scales = "free")

# Posterior ---------------------------------------------------------------
posterior_samples.ind = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  df.trial = df.behav %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-independent-data.wppl"),
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
save_data(posterior_samples.ind, paste(target_dir, "posterior_samples.rds", sep=FS))

# Chain Plots
df <- posterior_samples.ind %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_longer(cols=c(-Iteration, -Chain, -id), names_to = "param")
plots_chains = map(c("alpha", "gamma", "shape1", "shape2"), function(par){
  plots = map(c("blue", "green", "delta"), function(p) {
    fn <- paste(par, p, sep="_")
    p <- df %>% 
      filter(startsWith(param, fn)) %>% 
      ggplot(aes(x=Iteration, y = value, color=Chain, group = Chain)) + 
      geom_line() + 
      facet_wrap(id~param, scales="free", 
                 labeller = labeller(id = trial_names))
    
    ggsave(paste(target_dir, FS, paste("chains_", fn, ".png", sep = ""), sep=""), p)
    
    return(p)                 
  })
  return(plots)
})

# Mean posterior values for likelihood function
pars <- c("alpha_blue", "gamma_blue", "shape1_blue", "shape2_blue", 
          "alpha_green", "gamma_green", "shape1_green", "shape2_green",
          "alpha_delta", "gamma_delta", "shape1_delta", "shape2_delta", "theta_p")
evs.posterior.ind = posterior_samples.ind %>% 
  pivot_longer(cols = all_of(pars), names_to = "Parameter", values_to = "value") %>% 
  group_by(id, Parameter) %>%
  summarize(ev = mean(value), .groups = "drop_last") %>% 
  pivot_wider(names_from = "Parameter", values_from = "ev") %>% 
  # (order in thesis)
  dplyr::select(id, alpha_delta, gamma_delta, shape1_delta, shape2_delta,
                alpha_blue, gamma_blue, shape1_blue, shape2_blue, 
                alpha_green, gamma_green, shape1_green, shape2_green, 
                theta_p)

save_data(evs.posterior.ind, 
          paste(target_dir, "evs-posterior-independent-data.rds", sep=FS))


# Generate new independent tables -----------------------------------------
# generate set of dependent tables for each sample from posterior distribution (parameters)
sampled_tables = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  samples_trial <- posterior_samples.ind %>% filter(id == trial_id) %>% 
    filter(Iteration <= 100) #just use the first 100 samples, not all
  #samples_trial <- evs.posterior.ind %>% filter(id == trial_id) #use evs
  data_trial <- df.behav %>% filter(id == trial_id)
  map(seq(1, nrow(samples_trial)), function(idx){
    posterior_params <- samples_trial[idx,]
    data_webppl <- list(probs = data_trial, evs_params = posterior_params)
    
    samples.ind_tables <- webppl(
      program_file = here("webppl-model", "posterior-independent-data.wppl"),
      random_seed = params$seed_webppl, #?
      data_var = "data",
      model_var = "sample_table",
      data = data_webppl,
      packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
      inference_opts = list(method = "forward", samples = nrow(data_trial))
    ) %>% as_tibble() %>% 
      pivot_wider(names_from = "Parameter", values_from = "value") %>% 
      add_column(id = trial_id)
  })
})
save_data(sampled_tables, 
          paste(target_dir, "sampled-tables-posterior-predictive.rds", sep=FS))

# plot new tables sampled for each independent context with observed data
pp_plots_new_tables <- map(ind_trials, function(trial_id){
  df.trial <- df.ind %>% filter(id == trial_id)
  tit <- trial_names[[trial_id]]
  fn <- paste(target_dir, FS, paste("pp-tables-evs-posterior-", trial_id, 
                                    ".png", sep=""), sep="")
  p <- plot_new_tables(df.trial, sampled_tables, tit, fn)
})

# Posterior predictive ----------------------------------------------------
# Log likelihood plots
likelihood_fn = function(df.samples, df.grp){
  message(paste(df.grp$id))
  data_webppl <- list(probs = df.behav %>% filter(id == df.grp$id),
                      samples_posterior = df.samples)
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-independent-data.wppl"),
    random_seed = params$seed_webppl,
    data_var = "data",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS))
  ) 
  
  result = tibble(
    x_blue = samples.pp$x_blue,
    x_green = samples.pp$x_green,
    x_bg = samples.pp$x_bg,
    ll_X_new = samples.pp$ll_X, 
    ll_obs = samples.pp$ll_obs, 
    ll_X_blue = samples.pp$ll_x_blue, 
    ll_X_green = samples.pp$ll_x_green, 
    ll_X_bg = samples.pp$ll_x_bg, 
    ll_obs_blue = samples.pp$ll_obs_blue, 
    ll_obs_green = samples.pp$ll_obs_green, 
    ll_obs_bg = samples.pp$ll_obs_bg
  )
  return(result %>% add_column(id = df.grp$id))
}
pp_samples_ll = group_map(posterior_samples.ind %>% group_by(id),
                          likelihood_fn) %>% bind_rows()

# overall log likelihood of data from posterior predictive
ll_X_obs.mean = pp_samples_ll %>% group_by(id) %>% 
  summarize(ev = mean(ll_obs), .groups = "keep")

p <- pp_samples_ll %>% 
  ggplot(aes(x = ll_X_new)) + geom_density() +
  facet_wrap(~id, scales = "free", labeller = labeller(id = trial_names)) +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p
ggsave(paste(target_dir, "pp-log-likelihood-independent.png", sep = FS), p)

# log likelihood seperate for P(b), P(g) and P(b,g)
pp_samples_ll.long <- pp_samples_ll %>%
  pivot_longer(cols = c(-id,-x_blue, -x_green, -x_bg), 
               names_to = "prob", names_prefix = "ll_", values_to = "ll")
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
  facet_wrap(prob~id, ncol=5, scales = "free", 
             labeller = labeller(id = trial_names)) +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p
ggsave(paste(target_dir, "pp-log-likelihood-independent-separate.png", sep = FS), p)

