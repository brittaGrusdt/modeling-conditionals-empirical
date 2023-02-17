library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggpubr)

#' @param active_config str 'cf_prior_match_tables' or 'cf_prior_match_kl'
generate_model_states = function(active_config){
  Sys.setenv(R_CONFIG_ACTIVE = active_config)
  params <- config::get()
  states = generate_tbls_default_prior(params)
  return(states)
}

# Default prior ----------------------------------------------------
generate_tbls_default_prior = function(params){
  cns = create_causal_nets(params$p_noise, params$p_cp, params$p_ant, 
                           params$rels_dep, params$p_a, params$p_c)
  params$causal_nets_dep = cns$dep
  params$causal_nets_ind = cns$ind
  
  df.default_prior <- webppl(
    program_file = here("webppl-model", "run-state-prior.wppl"),
    data_var = "data",
    data = params,
    packages = params$packages[1:2],
    random_seed = params$seed_webppl
  ) 
  states.cf_prior = df.default_prior$prior$support %>% as_tibble()
  return(states.cf_prior)
}


  
