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
library(bayesplot)
library(ggthemes)
library(tidybayes)
library(emmeans)
library(brms)

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

theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
# Data --------------------------------------------------------------------
active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "dependent-contexts", sep=FS)
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)

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
  rename(blue = `blue falls`, green = `green falls`)

# Dependent Trials --------------------------------------------------------
df.dep = data.pe %>%
  rename(if_bg = `if blue falls green falls`, 
         if_nbg = `if blue does not fall green falls`, 
         if_gb = `if green falls blue falls`, 
         if_ngb = `if green does not fall blue falls`) %>% 
  dplyr::select(prolific_id, id, if_bg, if_nbg, if_gb, if_ngb, blue, green, 
                AC, `A-C`, `-AC`, `-A-C`) %>% 
  filter(!str_detect(id, "independent")) %>% 
  mutate(relation = case_when(str_detect(id, "if2") ~ "if2",
                              str_detect(id, "if1") ~ "if1"))

dep_trials = c(df.dep$id %>% unique())

##########################################################
####          Dependent Trials                        ####   
#### fit Zero-inflated Beta P(a), P(c|a), P(c|¬a)     ####
##########################################################

posterior_samples.dep = map_dfr(dep_trials, function(trial_id){
  message(trial_id)
  df.trial = df.dep %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-dependent-data.wppl"),
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
          "alpha_if_nbg", "gamma_if_nbg", "shape1_if_nbg", "shape2_if_nbg",
          "alpha_if_bg", "gamma_if_bg", "shape1_if_bg", "shape2_if_bg")
evs.posterior.dep = posterior_samples.dep %>% 
  pivot_longer(cols = all_of(pars), names_to = "Parameter", values_to = "value") %>% 
  group_by(id, Parameter) %>%
  summarize(ev = mean(value), .groups = "drop_last") %>% 
  pivot_wider(names_from = "Parameter", values_from = "ev") %>% 
  # (order in thesis)
  dplyr::select(id,
                alpha_if_bg, gamma_if_bg, shape1_if_bg, shape2_if_bg,
                alpha_if_nbg, gamma_if_nbg, shape1_if_nbg, shape2_if_nbg,
                alpha_blue, gamma_blue, shape1_blue, shape2_blue)

save_data(evs.posterior.dep, 
          paste(target_dir, "evs-posterior-dependent-data.rds", sep=FS))

# Generate new dependent tables -------------------------------------------
# get Posterior predictive data
# generate set of dependent tables with expected values of posterior distributions
sampled_tables = map_dfr(dep_trials, function(trial_id){
  message(trial_id)
  
  samples_trial <- posterior_samples %>% fitler(id == trial_id)
  map(seq(1, nrow(samples_trial)), function(idx){
    
    posterior_params <- samples_trial[idx,]
    data_trial <- df.dep %>% filter(id == trial_id)
    data_webppl <- list(probs = data_trial,
                        #probs = df.dep %>% filter(id == trial_id),
                        evs_params = posterior_params)
                        #evs_params = evs.posterior.dep %>% filter(id == trial_id))
    
    samples.dep_tables <- webppl(
      program_file = here("webppl-model", "posterior-dependent-data.wppl"),
      random_seed = params$seed_webppl,
      data_var = "data",
      model_var = "sample_table",
      data = data_webppl,
      packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
      inference_opts = list(method = "forward", samples = nrow(data_trial))
                            #samples = 1000)
    ) %>% as_tibble() %>% 
      pivot_wider(names_from = "Parameter", values_from = "value") %>% 
      add_column(id = trial_id)
  })
})

# plot new tables sampled for each dependent context with observed data
df.samples <- sampled_tables %>% 
  pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), 
               names_to = "cell", values_to = "value")

df.behav.long = df.dep %>% 
  dplyr::select(AC, `A-C`, `-AC`, `-A-C`, prolific_id, id) %>% 
  pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), 
               names_to = "cell", values_to = "value") 

colors <- c("empirical" = "darkgreen", "sampled" = "firebrick")
plots_new_tables <- map(dep_trials, function(trial_id){
  
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
ll_fn <- function(df.samples, df.grp){
  message(paste(df.grp$id))
  data_webppl <- list(probs = df.dep %>% filter(id == df.grp$id),
                      samples_posterior = df.samples)
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-dependent-data.wppl"),
    random_seed = params$seed_webppl,
    data_var = "data",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS))
  ) 
  
  result = tibble(ll_X_new = samples.pp$ll_X, 
                  ll_obs = samples.pp$ll_obs, 
                  ll_X_blue = samples.pp$ll_x_blue, 
                  ll_X_if_bg = samples.pp$ll_x_if_bg, 
                  ll_X_if_nbg = samples.pp$ll_x_if_nbg, 
                  ll_obs_blue = samples.pp$ll_obs_blue, 
                  ll_obs_if_bg = samples.pp$ll_obs_if_bg, 
                  ll_obs_if_nbg = samples.pp$ll_obs_if_nbg)
  return(result %>% add_column(id = df.grp$id))
}
pp_samples_ll = group_map(posterior_samples.dep %>% group_by(id), ll_fn) %>%
  bind_rows()

# overall log likelihood of data from posterior predictive
id_names <- c("if1_hh" = "if1:HI", 
              "if1_lh" = "if1:HI", 
              "if1_u-Lh" = "if1:LI", 
              "if1_uh" = "if1:UI", 
              "if2_hl" = "if2:UL", 
              "if2_ll" = "if2:UL", 
              "if2_u-Ll" = "if2:UL", 
              "if2_ul" = "if2:UL"
)

ll_X_obs.mean.all = pp_samples_ll %>% group_by(id) %>% 
  summarize(ev = mean(ll_obs), .groups = "keep")
p.all <- pp_samples_ll %>% 
  ggplot(aes(x = ll_X_new)) + geom_density() +
  facet_wrap(~id, scales = "free", ncol = 4, labeller = labeller(id = id_names)) +
  geom_point(data = ll_X_obs.mean.all, aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 20)) +
  labs(x = "log likelihood", y = "density")
p.all

ggsave(paste(target_dir, "pp-log-likelihood-dependent.png", sep = FS), p.all)

# log likelihood seperate for P(b), P(g|b) and P(g|¬b)
pp_samples_ll.long <- pp_samples_ll %>%
  pivot_longer(cols = c(-id), names_to = "prob", names_prefix = "ll_", 
               values_to = "ll")
# ll pp-samples
df.ll_X = pp_samples_ll.long %>% 
  filter(startsWith(prob, "X_") & prob != "X_new") %>% 
  mutate(prob = str_replace(prob, "X_", ""))

# ev ll observed data given params from posterior, seperately for blue/if_bg/if_nbg
ll_X_obs.mean = pp_samples_ll.long %>% group_by(id, prob) %>% 
  filter(startsWith(prob, "obs_")) %>% 
  summarize(ev = mean(ll), .groups = "drop_last") %>% 
  mutate(prob = str_replace(prob, "obs_", ""))

p.if1 <- df.ll_X  %>% filter(startsWith(id, "if1")) %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(id~prob, ncol=3, scales = "free") +
  geom_point(data = ll_X_obs.mean %>% filter(startsWith(id, "if1")), 
             aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "log likelihood", y = "density")
p.if1


p.if2 <- df.ll_X %>% filter(startsWith(id, "if2")) %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(id~prob, ncol=3, scales = "free") +
  geom_point(data = ll_X_obs.mean %>% filter(startsWith(id, "if2")), 
             aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "log likelihood", y = "density")
p.if2



