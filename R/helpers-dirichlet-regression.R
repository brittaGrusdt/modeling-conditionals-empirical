#' df.predictors long format with linear predictor in column 'eta' and
#' predicted category in column 'predictor' and column 'rowid' 
#' (one rowid for each table, i.e. predicted_cats: bg, b, g, none)
#' @param ref_response_cat str
transform_to_probs = function(df.predictors, ref_response_cat){
  df <- df.predictors %>% 
    # pivot_longer(cols=c(-category, -rowid), 
    #              names_to = "predicted_cat", values_to="eta") %>% 
    mutate(p_eta = exp(eta)) %>% 
    group_by(rowid, predictor) %>% 
    mutate(denominator = 1 + sum(p_eta),
           p_refcat = 1 / denominator, 
           p_hat = p_refcat * p_eta) %>% 
    dplyr::select(-eta, -p_eta, -denominator)
  
  df.all_categories = bind_rows(
    df %>% dplyr::select(-p_hat) %>% rename(p_hat = p_refcat) %>% 
      mutate(response_cat = ref_response_cat) %>% 
      distinct_at(vars(c(rowid, response_cat, predictor)), .keep_all = T), 
    df %>% dplyr::select(-p_refcat)
  )
  return(df.all_categories)           
}



# Dependent contexts ------------------------------------------------------
# reference categories: pblue: high, relation: if1
get_linear_predictors_dep = function(predicted_cat, model_params){
  mu = paste('mu', predicted_cat, sep='')
  pars = model_params[startsWith(names(model_params), mu)] 
  df = tibble(
    response_cat = predicted_cat,
    if1_high = pars[[paste(mu, '_Intercept', sep='')]], # reference categories
    if1_low = if1_high + pars[[paste(mu, '_pbluelow', sep = '')]],
    if1_unc = if1_high + pars[[paste(mu, '_pblueunc', sep = '')]],
    if1_uncl = if1_high + pars[[paste(mu, '_pblueuncl', sep = '')]],
    
    if2_high = if1_high + pars[[paste(mu, '_relationif2', sep='')]], 
    if2_low = if1_high + pars[[paste(mu, '_pbluelow', sep = '')]] + 
      pars[[paste(mu, '_relationif2', sep = '')]] + 
      pars[[paste(mu, '_pbluelow:relationif2', sep = '')]],
    
    if2_unc = if1_high + pars[[paste(mu, '_pblueunc', sep = '')]] +
      pars[[paste(mu, '_relationif2', sep = '')]] + 
      pars[[paste(mu, '_pblueunc:relationif2', sep = '')]],
    
    if2_uncl = if1_high + pars[[paste(mu, '_pblueuncl', sep = '')]] +
      pars[[paste(mu, '_relationif2', sep = '')]] + 
      pars[[paste(mu, '_pblueuncl:relationif2', sep = '')]],
  )
  return(df)
}

# get linear predictors for all posterior draws
get_linear_predictors_samples_dep = function(draws, response_cat){
  mu = paste('b_mu', response_cat, sep='')
  df = draws %>% 
    transmute(.chain, .iteration,
      if1_high = get(paste(mu, '_Intercept', sep='')),
      if1_low = if1_high + get(paste(mu, '_pbluelow', sep = '')),
      if1_unc = if1_high + get(paste(mu, '_pblueunc', sep = '')),
      if1_uncl = if1_high + get(paste(mu, '_pblueuncl', sep = '')),
      
      if2_high = if1_high + get(paste(mu, '_relationif2', sep='')), 
      if2_low = if1_high + get(paste(mu, '_pbluelow', sep = '')) + 
        get(paste(mu, '_relationif2', sep = '')) + 
        get(paste(mu, '_pbluelow:relationif2', sep = '')),
      
      if2_unc = if1_high + get(paste(mu, '_pblueunc', sep = '')) +
        get(paste(mu, '_relationif2', sep = '')) + 
        get(paste(mu, '_pblueunc:relationif2', sep = '')),
      
      if2_uncl = if1_high + get(paste(mu, '_pblueuncl', sep = '')) +
        get(paste(mu, '_relationif2', sep = '')) + 
        get(paste(mu, '_pblueuncl:relationif2', sep = ''))
    ) %>% 
    add_column(response_cat = response_cat) %>% 
    rowid_to_column()
  return(df)
}


# Independent contexts ----------------------------------------------------
# reference categories: pblue: high, pgreen: high
get_linear_predictors_ind = function(predicted_cat, model_params){
  mu = paste('mu', predicted_cat, sep='')
  pars = model_params[startsWith(names(model_params), mu)] 
  df = tibble(
    response_cat = predicted_cat,
    ind_hh = pars[[paste(mu, '_Intercept', sep='')]], # reference categories
    ind_hl = ind_hh + pars[[paste(mu, '_pgreenlow', sep = '')]],
    ind_uh = ind_hh + pars[[paste(mu, '_pblueunc', sep = '')]],
    ind_ul = ind_hh + pars[[paste(mu, '_pblueunc', sep = '')]] + 
      pars[[paste(mu, '_pgreenlow', sep = '')]] +
      pars[[paste(mu, '_pblueunc:pgreenlow', sep = '')]],
    ind_ll = ind_hh + pars[[paste(mu, '_pbluelow', sep = '')]] + 
      pars[[paste(mu, '_pgreenlow', sep = '')]] +
      pars[[paste(mu, '_pbluelow:pgreenlow', sep = '')]]
  )
  return(df)
}

# get linear predictors for all posterior draws
# reference category: pgreen: high, pblue:high
get_linear_predictors_samples_ind = function(draws, response_cat){
  mu = paste('b_mu', response_cat, sep='')
  df = draws %>% 
    transmute(
      ind_hh = get(paste(mu, '_Intercept', sep='')),
      ind_hl = ind_hh + get(paste(mu, '_pgreenlow', sep = '')),
      ind_uh = ind_hh + get(paste(mu, '_pblueunc', sep = '')),
      ind_ul = ind_hh + get(paste(mu, '_pblueunc', sep = '')) + 
        get(paste(mu, '_pgreenlow', sep = '')) +
        get(paste(mu, '_pblueunc:pgreenlow', sep = '')),
      ind_ll = ind_hh + get(paste(mu, '_pbluelow', sep = '')) + 
        get(paste(mu, '_pgreenlow', sep = '')) +
        get(paste(mu, '_pbluelow:pgreenlow', sep = ''))
    ) %>% 
    add_column(response_cat = response_cat) %>% 
    rowid_to_column()
  return(df)
}






