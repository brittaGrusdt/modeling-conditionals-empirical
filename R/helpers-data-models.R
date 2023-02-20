plot_new_tables <- function(df.trial, samples.tbls, tit, target_path){
  N = nrow(df.trial)
  grp <- rep(factor(c("bg", "b¬g", "¬bg", "¬b¬g"),
                    levels = c("bg", "b¬g", "¬bg", "¬b¬g")), N)
  trial_behav.mat <- df.trial %>% dplyr::select(AC, `A-C`, `-AC`, `-A-C`) %>%
    as.matrix()
  trial_behav <- c(t(trial_behav.mat))
  
  trial_samples.mat <- samples.tbls %>% 
    dplyr::select(AC, `A-C`, `-AC`, `-A-C`) %>% as.matrix()
  N_rep <- dim(trial_samples.mat)[1]/N
  yrep <- matrix(data=NA, nrow=N_rep, N*4)

  # N_rep x N samples
  for(i in seq(0, N_rep-1)){
    i_low <- i*N + 1
    i_up <- (i+1)*N

    vec <- c(t(trial_samples.mat[i_low:i_up, ]))
    for(j in seq(1, length(vec))){
      yrep[i + 1, j] <- vec[j]
    }
  }
  p <- ppc_dens_overlay_grouped(y=trial_behav, yrep = yrep, group = grp) +
    theme(panel.spacing = unit(2, "lines")) + ggtitle(tit)
  
  ggsave(target_path, p) 
  return(p)
}

# iterate over trials to make posterior predictive plots
make_pp_plots_new_dependent_tables = function(df.dep, sampled_tables, trials, 
                                              target_dir, fn_prefix){
  pp_plots_new_tables <- map(trials, function(trial_id){
    df.trial <- df.dep %>% filter(id == trial_id)
    tit <- switch(trial_id,
                  "if1_hh"=expression("if"[1]*":HI"),
                  "if1_uh"=expression(paste(`if`[1], ":UI")),
                  "if1_u-Lh"=expression("if"[1]*":U"^-{}*"I"),
                  "if1_lh"=expression("if"[1]*":LI"),
                  "if2_hl"=expression("if"[2]*":HL"),
                  "if2_ul"=expression("if"[2]*":UL"),
                  "if2_u-Ll"=expression("if"[2]*":U"^-{}*"L"),
                  "if2_ll"=expression("if"[2]*":LL"))
    fn <- paste(fn_prefix, "_", trial_id, ".png", sep="")
    target_path <- paste(target_dir, FS, fn, sep="")
    df.samples <- sampled_tables %>% filter(id == trial_id)
    p <- plot_new_tables(df.trial, df.samples, tit, target_path)
    return(p)
  })
  return(pp_plots_new_tables)
}


sample_tables = function(df.behav, trials, params, fn_wppl_program, repetitions=1){
  all_samples <- map(seq(1, repetitions), function(i_rep) {
    sampled_tables = map_dfr(trials, function(trial_id){
      message(trial_id)
      samples_trial <- params %>% filter(id == trial_id) 
      data_trial <- df.behav %>% filter(id == trial_id)
      
      map(seq(1, nrow(samples_trial)), function(idx){
        posterior_params <- samples_trial[idx,]
        data_webppl <- list(probs = data_trial,
                            evs_params = posterior_params)
        
        samples.tbls <- webppl(
          program_file = fn_wppl_program, #here("webppl-model", "posterior-dependent-data.wppl"),
          data_var = "data",
          model_var = "sample_table",
          data = data_webppl,
          packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
          inference_opts = list(method = "forward", samples = nrow(data_trial))
        ) %>% as_tibble() %>% 
          pivot_wider(names_from = "Parameter", values_from = "value") %>% 
          add_column(id = trial_id)
        return(samples.tbls)
      })
    })
    return(sampled_tables)
  }) %>% bind_rows()
  return(all_samples)
}

simulate_zoib_data = function(N, alpha, gamma, shape1, shape2){
  zero_ones = rbernoulli(N, alpha)
  simulations = map_dbl(zero_ones, function(c){
    if(c) { p = rbernoulli(1, gamma)
    } else {
      p = rbeta(1, shape1, shape2)
    }
    return(p)
  })
  return(simulations)
}