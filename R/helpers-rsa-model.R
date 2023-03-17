get_assertable_utterances = function(tbls.wide, theta){
  tbls.wide %>% 
    rename(bg = AC, b = `A-C`, g = `-AC`, none = `-A-C`) %>% 
    add_probs() %>% 
    mutate(p_AC = bg, `p_A-C` = b, `p_-AC` = g, `p_-A-C` = none) %>% 
    pivot_longer(starts_with("p_"), names_to = "utt", values_to = "prob") %>% 
    mutate(utt = case_when(utt == "p_AC" ~ standardized.sentences$bg, 
                           utt == "p_A-C" ~ standardized.sentences$b, 
                           utt == "p_-AC" ~ standardized.sentences$g,
                           utt == "p_-A-C" ~ standardized.sentences$none, 
                           utt == "p_a" ~ standardized.sentences$only_b,
                           utt == "p_c" ~ standardized.sentences$only_g,
                           utt == "p_na" ~ standardized.sentences$only_nb,
                           utt == "p_nc" ~ standardized.sentences$only_ng, 
                           utt == "p_c_given_a" ~ standardized.sentences$if_bg,
                           utt == "p_c_given_na" ~ standardized.sentences$if_nbg,
                           utt == "p_nc_given_a" ~ standardized.sentences$if_bng,
                           utt == "p_nc_given_na" ~ standardized.sentences$if_nbng,
                           utt == "p_a_given_c" ~ standardized.sentences$if_gb,
                           utt == "p_a_given_nc" ~ standardized.sentences$if_ngb,
                           utt == "p_na_given_c" ~ standardized.sentences$if_gnb,
                           utt == "p_na_given_nc" ~ standardized.sentences$if_ngnb, 
                           utt == "p_likely_a" ~ standardized.sentences$might_b,
                           utt == "p_likely_c" ~ standardized.sentences$might_g,
                           utt == "p_likely_na" ~ standardized.sentences$might_nb,
                           utt == "p_likely_nc" ~ standardized.sentences$might_ng
    )) %>% 
    mutate(u_assertable = case_when(str_detect(utt, "might") ~ T, 
                                    prob >= theta ~ T, T ~ F)) %>% 
    filter(u_assertable) %>% dplyr::select(-u_assertable)
}


# utt_costs is a tibble with columns 'utterance', 'cost'
prepare_data_for_wppl <- function(config_cns = "default_cns", 
                                  config_weights_relations = "semi_informative",
                                  config_fits = NA,
                                  utt_costs = NA,
                                  extra_packages = NA
                                  ){
  
  # retrieve default params
  Sys.setenv(R_CONFIG_ACTIVE = "pathes")
  params <- config::get()
  if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive=T)
  
  # add custom parameters as given by function arguments
  # 1. causal networks
  Sys.setenv(R_CONFIG_ACTIVE = config_cns)
  pars.cns <- config::get()
  cns = create_causal_nets(pars.cns$dep_noise, 
                           pars.cns$dep_causal_power, 
                           pars.cns$dep_marginal, 
                           pars.cns$rels_dep,
                           pars.cns$ind_a, 
                           pars.cns$ind_c)
  params$causal_nets_dep = cns$dep
  params$causal_nets_ind = cns$ind

  data.behav <- read_csv(here(params$dir_data, params$fn_cleaned_data)) %>% 
    dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
    translate_standardized2model() 
  
  # 2. add default rsa parameters and custom packages
  Sys.setenv(R_CONFIG_ACTIVE = "rsa_params")
  rsa_params <- config::get()
  params <- c(params, rsa_params)

  if(!is.na(extra_packages)){
    pathes_extra_packages <- map_chr(extra_packages, function(name_pck){
      return(params[[name_pck]])
    })
    params$packages <- c(params$packages, pathes_extra_packages)
  }
  
  # 3. retrieve observed data each participant and trial + ratios
  pars.observations <- get_observed_data(data.behav, params)
  # 4. draw samples from prior with cns specified above
  pars.rsa_states <- get_rsa_states(params)
  # 5. get parameters of data likelihood functions
  pars.likelihoods <- get_likelihood_params_fitted_data()
  
  # 6. join these into single list
  params <- c(params, pars.observations, pars.likelihoods, pars.rsa_states)
  
  # 7. custom utterance cost
  if(!is.na(utt_costs)){
    params$utt_cost <- utt_costs
  }
  
  # 8. prior relations (for context-weights)
  Sys.setenv(R_CONFIG_ACTIVE =  "priors_relations")
  par.relations <- config::get()
  params$prior_relations <- par.relations[[config_weights_relations]]
  
  if(!is.na(config_fits)){
    Sys.setenv(R_CONFIG_ACTIVE = config_fits)
    par = config::get()
    params <- c(params, par)
  }
  
  return(params)
}

# brings params in particular utterance cost into format expected by wppl
format_param_samples = function(samples){
  if(length(samples$Parameter[str_detect(samples$Parameter, "utt_cost")]) > 0) {
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
    ) 
  } else {
    samples.formatted <-  samples %>% 
      filter(!startsWith(as.character(Parameter), "utts.")) %>% 
      pivot_wider(names_from=Parameter, values_from=value)
  }
  df.res <- samples.formatted %>% 
    mutate(sample_id = paste("chain", Chain, "-Iteration", Iteration, sep="")) %>% 
    rowid_to_column()
  
  return(df.res)
}





