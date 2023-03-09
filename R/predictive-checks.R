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

path_model_file = paste(params$dir_wppl_code, "predictive-checks.wppl", sep=FS)

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
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()
params$prior_relations <- par_relations[["informative"]]
#params$utt_cost <- tibble(utterance = params$utterances, cost=params$utt_cost)



# Helper functions --------------------------------------------------------
format_param_samples = function(samples){
  samples.utt_cost <- samples %>% 
    filter(startsWith(as.character(Parameter), "utts.")) %>% 
    mutate(Parameter = as.character(Parameter), 
           Parameter = str_replace(Parameter, "utts.", "")) %>% 
    rename(cost = value, utterance = Parameter) %>% 
    group_by(Iteration, Chain) %>% 
    nest() %>% rename(utt_cost=data)
  
  samples.formatted <- left_join(
    samples.utt_cost,
    samples %>% 
      filter(!startsWith(as.character(Parameter), "utts.")) %>% 
      pivot_wider(names_from=Parameter, values_from=value)
  ) %>% 
    mutate(sample_id = paste("chain", Chain, "-Iteration", Iteration, sep=""))
  return(samples.formatted)
}

# Prior predictive --------------------------------------------------------
params$par_fit <- c("alpha", "theta", "utt_cost") #, "gamma")
model <- "var prior = function(){
  return(priorSample(data.par_fit))
}
"
# get samples from prior distribution
prior_samples <- webppl(program_code = model, 
                        data_var = "data",
                        model_var = "prior",
                        data = params,
                        inference_opts = list(method = "forward",
                                              samples = 4000),
                        packages = params$packages) %>% as_tibble() %>% 
  mutate(Parameter = as.character(Parameter), 
         Parameter = str_replace(Parameter, "utt_cost", "utts"))
params$sampled_params <- format_param_samples(prior_samples)[1:2,]

# then run RSA-model once with each sampled set of parameters
data <- webppl(program_file = path_model_file,
               data = params,
               data_var = "data",
               random_seed = params$seed_webppl,
               packages = params$packages)

prior_predictive <- data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows()


# Posterior predictive ----------------------------------------------------
fn <- "alpha_theta_utt_cost"
# get samples from posterior distribution
posterior_samples <- readRDS(here(params$dir_results, 
                             paste("mcmc-posterior-fit-context-predictions-",
                                   fn, ".rds", sep="")))
params$sampled_params <- format_param_samples(posterior_samples)[1:2,]

# then run RSA-model once with each sampled set of parameters
data <- webppl(program_file = path_model_file,
               data = params,
               data_var = "data",
               random_seed = params$seed_webppl,
               packages = params$packages)

posterior_predictive <- data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows()



