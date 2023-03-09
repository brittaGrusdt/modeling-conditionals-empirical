library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

# parameters to run webppl model
params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
result_dir = params$dir_results

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
#path_model_file = paste(params$dir_wppl_code, "model-single-run-by-contexts.wppl", sep=FS)
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()
params$prior_relations <- par_relations[["informative"]]
#params$utt_cost <- tibble(utterance = params$utterances, cost=params$utt_cost)


# Prior predictive --------------------------------------------------------


# Posterior predictive ----------------------------------------------------
fn <- "alpha_theta_utt_cost"
posterior_samples <- readRDS(here(params$dir_results, 
                             paste("mcmc-posterior-fit-context-predictions-",
                                   fn, ".rds", sep="")))
samples.utt_cost <- posterior_samples %>% 
  filter(startsWith(as.character(Parameter), "utts.")) %>% 
  mutate(Parameter = as.character(Parameter), 
         Parameter = str_replace(Parameter, "utts.", "")) %>% 
  rename(cost = value, utterance = Parameter) %>% 
  group_by(Iteration, Chain) %>% 
  nest() %>% rename(utt_cost=data)

samples <- left_join(
  posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  pivot_wider(names_from=Parameter, values_from=value), 
  samples.utt_cost
) %>% 
  mutate(sample_id = paste("chain", Chain, "-Iteration", Iteration,))

params$sampled_params <- samples[1:2,]
path_model_file = paste(params$dir_wppl_code, "predictive-checks.wppl", sep=FS)


data <- webppl(program_file = path_wppl_file,
               data = params,
               data_var = "data",
               random_seed = params$seed_webppl,
               packages = params$packages)

pp_samples <- data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows()



