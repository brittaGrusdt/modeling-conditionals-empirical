plot_new_tables <- function(df.trial, sampled_tables, tit, target_path){
  N = nrow(df.trial)
  grp <- rep(factor(c("bg", "b¬g", "¬bg", "¬b¬g"),
                    levels = c("bg", "b¬g", "¬bg", "¬b¬g")), N)
  trial_behav.mat <- df.trial %>% dplyr::select(AC, `A-C`, `-AC`, `-A-C`) %>%
    as.matrix()
  trial_behav <- c(t(trial_behav.mat))
  
  trial_samples.mat <- sampled_tables %>% 
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
