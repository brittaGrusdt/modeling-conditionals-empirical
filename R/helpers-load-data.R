# load likelihood parameters fitted to data
get_likelihood_params_fitted_data = function(params){
  pars.dep <- readRDS(here(params$dir_dep_likelihoods, params$fn_likelihoods)) %>% 
    dplyr::select(variable, mean, id) %>% 
    rename(p = variable) %>% 
    mutate(p = str_replace(p, "if_", "if")) %>% 
    separate(p, into = c("param", "p"), sep="_") %>% 
    mutate(p = str_replace(p, "if", "if_")) %>% 
    pivot_wider(names_from = "param", values_from = "mean")
  
  pars.ind <- readRDS(here(params$dir_ind_likelihoods, params$fn_likelihoods))
  pars.ind.zoib <- pars.ind  %>% 
    dplyr::select(variable, mean, id) %>%
    filter(variable != "sd_delta") %>% 
    rename(p = variable) %>% 
    separate(p, into = c("param", "p"), sep="_") %>% 
    pivot_wider(names_from = "param", values_from = "mean")
  
  pars.ind.gaussian <- pars.ind %>% 
    dplyr::select(variable, mean, id) %>% 
    rename(p = variable) %>% 
    filter(p == "sd_delta") %>% 
    separate(p, into = c("param", "p"), sep="_") %>% 
    pivot_wider(names_from = "param", values_from = "mean")
  
  pars = list(likelihoods_zoib = bind_rows(pars.dep, pars.ind.zoib),
              likelihoods_gaussian = pars.ind.gaussian)
  return(pars)
}


draw_states_from_prior = function(params){
  prior <- webppl(program_file = here("webppl-model", "run-state-prior.wppl"),
                  data_var = "data",
                  data = params,
                  random_seed = params$seed_webppl,
                  packages = c("webppl-model/node_modules/conditionalsHelpers", 
                               "webppl-model/node_modules/conditionalsDefault")) 
  states = prior$prior$support %>% as_tibble()
  return(states)
}

get_rsa_states = function(params){
  states <- draw_states_from_prior(params)
  pars <- list(prior_samples = states)
  return(pars)
}

get_observed_data = function(data.behav, params){
  data.uc = data.behav %>%  filter(!is.na(uc_task)) %>% 
    dplyr::select(prolific_id, id, utterance) %>% 
    add_column(n = 1) %>% group_by(id, utterance)
  data.pe = data.behav %>% 
    dplyr::select(prolific_id, id, utt.standardized, pe_task) %>% 
    pivot_wider(names_from = "utt.standardized", values_from = "pe_task") %>% 
    rename(AC = `both blocks fall`, 
           `A-C` = `blue falls but green does not fall`, 
           `-AC` = `green falls but blue does not fall`, 
           `-A-C` = `neither block falls`)
  # all subjects and trials
  data.observed = left_join(data.uc, data.pe)
  
  # observed ratios per context
  contexts <- data.behav$id %>% unique()
  ratios.observed <- left_join(
    tibble(id = contexts, utterance = list(params$utterances)) %>%
      unnest(c(utterance)),
    data.observed %>% group_by(id, utterance) %>% count(), 
  ) %>% mutate_if(is.numeric, coalesce, 0) %>% 
    group_by(id) %>% 
    mutate(ratio = n / sum(n))
  
  pars <- list(observations = data.observed, 
               observed_utts_ratios = ratios.observed)
  return(pars)
}

