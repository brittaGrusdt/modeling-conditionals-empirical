# load likelihood parameters fitted to data
get_likelihood_params_fitted_data = function(result_dir, params){
  
  dep_data_dir = paste(result_dir, "dependent-contexts", sep=FS)
  ind_data_dir = paste(result_dir, "independent-contexts", sep=FS)
  
  pars.dep <- readRDS(paste(dep_data_dir, "evs-posterior-dependent-data.rds", sep=FS)) %>% 
    pivot_longer(cols=c(-id), names_to = "p", values_to = "value") %>% 
    mutate(p = str_replace(p, "if_", "if")) %>% 
    separate(p, into = c("param", "p"), sep="_") %>% 
    mutate(p = str_replace(p, "if", "if_")) %>% 
    pivot_wider(names_from = "param", values_from = "value")
  
  pars.ind <- readRDS(paste(ind_data_dir, "evs-posterior-independent-data.rds", sep=FS)) %>% 
    dplyr::select(-theta_p) %>% 
    pivot_longer(cols=c(-id), names_to = "p", values_to = "value") %>% 
    separate(p, into = c("param", "p"), sep="_") %>% 
    pivot_wider(names_from = "param", values_from = "value")
  
  # noise added vs.subtracted to P(b)*P(g)
  pars.ind.theta_p <- readRDS(
    paste(ind_data_dir, "evs-posterior-independent-data.rds", sep=FS)
  ) %>% dplyr::select(id, theta_p) %>%
    rename(theta = theta_p)
  
  pars = list(likelihoods_zoib = bind_rows(pars.dep, pars.ind),
              likelihoods_bernoulli = pars.ind.theta_p)
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

