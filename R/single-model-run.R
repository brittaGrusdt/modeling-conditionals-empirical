library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "plot-functions.R"))
source(here("R", "helpers-load-data.R"))

theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
# use special causal nets, with different marginal distributions than uniform
# cns = create_causal_nets(params$p_noise, params$p_cp, params$p_ant,
#                          params$rels_dep, params$p_a, params$p_c)
# params$causal_nets_dep = cns$dep
# params$causal_nets_ind = cns$ind

# behavioral data
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
  translate_standardized2model() 

# add parameters to run webppl model
# observed data each participant and trial + ratios
pars.observations <- get_observed_data(data.behav, params)
# draw samples from prior
pars.rsa_states <- get_rsa_states(params)
# load likelihood parameters fitted to data
pars.likelihoods <- get_likelihood_params_fitted_data(params)
# if not provided, all utterances equally likely
# params$p_utts = rep(1 / length(params$utterances), length(params$utterances))
params <- c(params, pars.observations, pars.likelihoods, pars.rsa_states)

# Run Model ---------------------------------------------------------------
# 1. predictions by contexts
path_model_file = paste(params$dir_wppl_code, "model-single-run-by-contexts.wppl", sep=FS)
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()
params$prior_relations <- par_relations[["informative"]]

posterior <- run_webppl(path_model_file, params)
model.predictions <- posterior %>% 
  map(function(x){as_tibble(x) %>% mutate(ll_ci = as.numeric(ll_ci))}) %>% 
  bind_rows() %>% group_by(id) %>% 
  mutate(p_hat_round = round(p_hat, 2)) %>% 
  arrange(desc(p_hat))

production.joint = left_join(model.predictions, 
                             params$observed_utts_ratios %>% rename(p = ratio)) %>% 
  rename(model = p_hat, behavioral = p) %>% 
  mutate(relation = case_when(startsWith(id, "independent") ~ "independent", 
                              startsWith(id, "if1") ~ "if1",
                              startsWith(id, "if2") ~ "if2"))

plot_correlation(production.joint, color = "id") + facet_wrap(~id)














