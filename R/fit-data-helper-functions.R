
#'  Fit zero-one inflated beta distribution to data
#'  
#'  `evs.posterior` returns the expected value for each fitted parameter
#'  
#'  @param data data frame with column 'id' containing names of trials.
#'  @param trials character vector with trial names. 
#'  "all": data not filtered
#'  "if1": only data with relation=if1 used
#'  "if2": only data with relation=if2 used
#'  @param probs character vector whose elements are columns of 'data'
#'  @returns  `evs.posterior` 
#'  @examples tba
fit_zoib = function(data, trials, probs){
  
  samples.posterior = map_dfr(trials, function(trial_id){
    message(trial_id)
    if(trial_id == "all") {
      df.trial = data
    } else if(trial_id == "if1") {
      df.trial = data %>% filter(relation == "if1")
    } else if(trial_id == "if2") {
      df.trial = data %>% filter(relation == "if2")
    } else {
      df.trial = data %>% filter(id == trial_id)
    }
    
    samples.posterior_trial = map_dfr(probs, function(p) {
      df.trial_p <- df.trial %>% rename(p = !!p) %>% filter(!is.na(p))
      data_webppl = list(probs = df.trial_p)
      samples.posterior <- webppl(
        program_file = here("webppl-model", "posterior-zoib-data.wppl"),
        data_var = "data",
        model_var = "non_normalized_posterior",
        data = data_webppl,
        inference_opts = list(method = "MCMC",
                              samples = 1000,
                              lag = 2,
                              burn = 10000,
                              verbose = T),
        chains = 4, 
        cores = 4) %>% as_tibble() %>% 
        pivot_wider(names_from = "Parameter", values_from = "value") %>% 
        add_column(id = trial_id, p = p)
      
      # evs = samples.posterior %>% group_by(Parameter) %>% 
      #   summarize(ev = mean(value)) %>% 
      #   pivot_wider(names_from = "Parameter", values_from = "ev")
      
      # return(evs)
      return(samples.posterior)
    })
    return(samples.posterior_trial)
  })
  return(samples.posterior)
}


get_evs_zoib_samples = function(samples) {
  evs = samples %>% 
    pivot_longer(cols = c(shape1, shape2, gamma, alpha), 
                 names_to = "Parameter", 
                 values_to = "value") %>% 
    group_by(id, p, Parameter) %>%
    summarize(ev = mean(value), .groups = "drop_last") %>%
    pivot_wider(names_from = "Parameter", values_from = "ev")
  return(evs)
}


#' 
#' @param samples data frame with column p containing the fitted probabilities
#' @param data.behav PE-task data
posterior_predictive_zoib = function(samples, data.behav){
  
  pp.samples = group_map(samples %>% group_by(id, p), function(df.samples, df.grp){
    message(paste(df.grp$id, df.grp$p))
    # get empirical data for id-probability combination
    if(df.grp$id == "all") {
      df.trial = data.behav
    } else if(df.grp$id == "if1") {
      df.trial = data.behav %>% filter(relation == "if1")
    } else if(df.grp$id == "if2") {
      df.trial = data.behav %>% filter(relation == "if2")
    } else {
      df.trial = data.behav %>% filter(id == df.grp$id)
    }

    data_webppl = list(probs = df.trial[[df.grp$p]], 
                       samples_posterior = list(shape1 = df.samples$shape1,
                                                shape2 = df.samples$shape2,
                                                alpha = df.samples$alpha,
                                                gamma = df.samples$gamma))
    samples.pp <- webppl(
      program_file = here("webppl-model", "posterior-predictive-zoib-data.wppl"),
      data_var = "data",
      data = data_webppl
    )
    # byrow is important!! by default: columns are used
    X_new <- unlist(samples.pp$X_new) %>%
      matrix(ncol = length(df.trial[[df.grp$p]]), byrow=T)
    ll_X_new <- samples.pp$ll_X_new
    ll_X_obs = samples.pp$ll_X_obs
    
    result = tibble(X_new = list(X_new), 
                    ll_X_new = list(ll_X_new), 
                    ll_X_obs = list(ll_X_obs))
    return(result %>% add_column(id = df.grp$id, p = df.grp$p))
  }) %>% bind_rows()
  
  return(pp.samples)
} 


#'  returns density of log likelihoods for samples from posterior predictive
#'  with expected value of the log likelihood for observed data
#'
#' @param samples.pp data frame with columns id, p, X_new, ll_X_new, ll_X_obs
#' @param p_cols named vector from entries of column p to color-string
#' 
plot_pp_ll = function(samples.pp, p_cols){
  
  ll_data = group_map(samples.pp %>% group_by(id, p), function(df.samples, df.grp){
    
    n_samples = unlist(df.samples$ll_X_new) %>% length()
    pp_data = unlist(df.samples$X_new)
    n_col = pp_data %>% length() / n_samples # nb. of data points for this trial
    
    X_new <- pp_data %>% matrix(ncol = n_col, byrow=T)
    ll_X_new <- unlist(df.samples$ll_X_new)
    ll_X_obs = unlist(df.samples$ll_X_obs)
  
    df <- tibble(ll_X = ll_X_new, ll_X_obs = ll_X_obs) %>% 
      add_column(id = df.grp$id, p = df.grp$p) %>% 
      mutate(ll_X_obs = as.numeric(ll_X_obs))
    return(df)
  }) %>% bind_rows() %>% group_by(id, p)
  
  ll_X_obs.mean = ll_data %>% summarize(ev = mean(ll_X_obs), .groups = "keep")
  
  p <- ll_data %>% 
    ggplot(aes(x = ll_X, color = p)) + 
    geom_density() +
    facet_grid(p~id, scales = "free") +
    geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2) +
    scale_color_manual(name = "probability", values = p_cols) + 
    theme(legend.position = "top") +
    labs(x = "log likelihood for samples from posterior predictive with expected log likelihood of observed data", 
         y = "density")
  
  return(p)
}






