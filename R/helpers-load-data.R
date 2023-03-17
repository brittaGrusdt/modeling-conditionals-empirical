# load likelihood parameters fitted to data
get_likelihood_params_fitted_data = function(){
  Sys.setenv(R_CONFIG_ACTIVE = "pathes")
  params <- config::get()
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
  wppl_code <- "
    globalStore.n_prior_samples = data['nb_rsa_states'][0]
    var prior = state_prior()
    var distributions = {prior:  prior}
    distributions
  "
  Sys.setenv(R_CONFIG_ACTIVE = "rsa_params")
  params <- c(params, config::get())
  
  Sys.setenv(R_CONFIG_ACTIVE = "pathes")
  pathes <- config::get()
  packages <- c(pathes$conditionalsHelpers, pathes$conditionalsDefault)
  
  prior <- webppl(program_code = wppl_code,
                  data_var = "data",
                  data = params,
                  random_seed = params$seed_webppl,
                  packages = packages) 
  states = prior$prior$support %>% as_tibble()
  return(states)
}

get_rsa_states = function(params){
  states <- draw_states_from_prior(params)
  pars <- list(prior_samples = states)
  return(pars)
}

# if cutoff is TRUE (default=FALSE), utterances that were overall not selected 
# (when proportion rounded to two digits = 0) are filtered out
get_observed_data = function(data.behav, params, cutoff=F){
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
  if(cutoff){
    utts_low_freq = data.observed %>% 
      group_by(utterance) %>% dplyr::count() %>% 
      arrange(desc(n)) %>% ungroup() %>% mutate(N=sum(n), ratio=round(n/N, 2)) %>% 
      filter(ratio==0) %>% pull(utterance)
    data.observed <- data.observed %>% filter(!utterance %in% utts_low_freq)
  }
  
  # observed ratios per context
  contexts <- data.behav$id %>% unique()
  if(is.null(params$utterances)) stop(message("argument parameters needs entry 'utterances' containing all MODEL utterances"))
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

