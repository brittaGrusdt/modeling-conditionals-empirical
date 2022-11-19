library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggpubr)

#' @param active_config str one of 'context_sensitive_prior', 
#' 'abstract_default_params'
generate_model_states = function(active_config){
  Sys.setenv(R_CONFIG_ACTIVE = active_config)
  params <- config::get()
  if(active_config == "context_sensitive_prior") {
    states = generate_tbls_context_sensitive(params) 
  } else if(active_config == "abstract_default_params") {
    states = generate_tbls_abstract_prior(params)
  }
  return(states)
}

# Context-sensitive prior -------------------------------------------------
generate_tbls_context_sensitive = function(params){
  params.fitted_dirichlet = read_csv(
    here(params$dir_data, params$fn_fitted_alphas_by_context)
  )
  states = map_dfr(seq(1, nrow(params.fitted_dirichlet)), function(row){
    print(params.fitted_dirichlet[row, ]$id)
    params$alphas <- params.fitted_dirichlet[row,2:5] %>% as.numeric()
    df.context_sensitive_prior <- webppl(
      program_file = here("webppl-model", "run-state-prior.wppl"),
      data_var = "data",
      data = params,
      packages = params$packages[1])
    df.context_sensitive_prior$prior$support %>% as_tibble() %>% 
      rowid_to_column("bn_id") %>% group_by(bn_id) %>% 
      unnest(c(probs, support)) %>% 
      pivot_wider(names_from = "support", values_from = "probs") %>% 
      add_column(id = params.fitted_dirichlet[row, ]$id) %>% 
      mutate(trial = id, tmp = paste(AC, `A-C`, `-AC`, `-A-C`, collapse = "_")) %>% 
      unite("bn_id", trial, tmp, sep = "_")
  }) %>% 
  mutate(r = id, cn = "", 
         table.support=list(c("AC", "A-C", "-AC", "-A-C")),
         table.probs = list(c(AC, `A-C`, `-AC`, `-A-C`))) %>% 
    dplyr::select(-AC, -`A-C`, -`-AC`, -`-A-C`)
  return(states)
}
# Abstract state prior ----------------------------------------------------
generate_tbls_abstract_prior = function(params){
  cns = create_causal_nets(params$p_noise, params$p_cp, params$p_ant, 
                           params$rels_dep, params$p_a, params$p_c)
  params$causal_nets_dep = cns$dep
  params$causal_nets_ind = cns$ind
  
  df.abstract_state_prior <- webppl(
    program_file = here("webppl-model", "run-state-prior.wppl"),
    data_var = "data",
    data = params,
    packages = params$packages[1:2]) 
  states.abstract_prior = df.abstract_state_prior$prior$support %>%
    as_tibble()
  return(states.abstract_prior)
}


  
