library(here)
library(tidyverse)
library(tibble)
library(broom)
library(broom.mixed)
library(knitr)
library(ExpDataWrangling)
library(ModelUtils)
library(rwebppl)
library(bayesplot)
library(xtable)
library(ggpubr)
library(tidyselect)
library(bayesplot)
library(ggthemes)
library(tidybayes)
library(emmeans)
library(brms)

source(here("R", "helpers-data-models.R"))

# for plots
theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
trial_names <- c("if1_hh"=expression("if"[1]*":HI"),
                 "if1_uh"=expression(paste(`if`[1], ":UI")),
                 "if1_u-Lh"=expression("if"[1]*":U"^-{}*"I"),
                 "if1_lh"=expression("if"[1]*":LI"),
                 "if2_hl"=expression("if"[2]*":HL"),
                 "if2_ul"=expression("if"[2]*":UL"),
                 "if2_u-Ll"=expression("if"[2]*":U"^-{}*"L"),
                 "if2_ll"=expression("if"[2]*":LL"))
prob_names <- c("blue"="P(b)", "if_bg" = "P(g|b)", "if_nbg"="P(g|¬b)")

# Data --------------------------------------------------------------------
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "dependent-contexts",
                   "zoib-model", sep=FS)
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)

data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
  translate_standardized2model() 

data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task) %>% 
  pivot_wider(names_from = "utt.standardized", values_from = "pe_task") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`) %>% 
  get_controlled_factors() %>% dplyr::select(-relation_type) %>% 
  rename(blue = `blue falls`, green = `green falls`)

# Dependent Trials --------------------------------------------------------
df.dep = data.pe %>%
  rename(if_bg = `if blue falls green falls`, 
         if_nbg = `if blue does not fall green falls`, 
         if_gb = `if green falls blue falls`, 
         if_ngb = `if green does not fall blue falls`) %>% 
  dplyr::select(prolific_id, id, if_bg, if_nbg, if_gb, if_ngb, blue, green, 
                AC, `A-C`, `-AC`, `-A-C`) %>% 
  filter(!str_detect(id, "independent")) %>% 
  mutate(relation = case_when(str_detect(id, "if2") ~ "if2",
                              str_detect(id, "if1") ~ "if1"))

dep_trials = c(df.dep$id %>% unique())

##########################################################
####          Dependent Trials                        ####   
#### fit Zero-inflated Beta P(a), P(c|a), P(c|¬a)     ####
##########################################################

posterior_samples.dep = map_dfr(dep_trials, function(trial_id){
  message(trial_id)
  df.trial = df.dep %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  
  samples.posterior <- webppl(
    program_file = here("webppl-model", "posterior-dependent-data.wppl"),
    data_var = "data",
    model_var = "non_normalized_posterior",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
    inference_opts = list(method = "MCMC",
                          samples = 1000,
                          lag = 2,
                          burn = 10000,
                          verbose = T),
    chains = 4, 
    cores = 4) %>%
    as_tibble() %>% 
    pivot_wider(names_from = "Parameter", values_from = "value") %>% 
    add_column(id = trial_id)
  
  return(samples.posterior)
})
save_data(posterior_samples.dep, paste(target_dir, "posterior_samples.rds", sep=FS))

# Chain Plots
df <- posterior_samples.dep %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_longer(cols=c(-Iteration, -Chain, -id), names_to = "param")

plots_chains = map(c("alpha", "gamma", "shape1", "shape2"), function(par){
  plots = map(c("blue", "if_bg", "if_nbg"), function(p) {
    fn <- paste(par, p, sep="_")
    p <- df %>% 
      filter(startsWith(param, fn)) %>% 
      ggplot(aes(x=Iteration, y = value, color=Chain, group = Chain)) + 
      geom_line() + facet_wrap(id~param, scales="free")
    
    ggsave(paste(target_dir, FS, paste("chains_", fn, ".png", sep = ""), sep=""), p)
    
    return(p)                 
  })
  return(plots)
})

pars <- c("alpha_blue", "gamma_blue", "shape1_blue", "shape2_blue", 
          "alpha_if_nbg", "gamma_if_nbg", "shape1_if_nbg", "shape2_if_nbg",
          "alpha_if_bg", "gamma_if_bg", "shape1_if_bg", "shape2_if_bg")
evs.posterior.dep = posterior_samples.dep %>% 
  pivot_longer(cols = all_of(pars), names_to = "Parameter", values_to = "value") %>% 
  group_by(id, Parameter) %>%
  summarize(ev = mean(value), .groups = "drop_last") %>% 
  pivot_wider(names_from = "Parameter", values_from = "ev") %>% 
  # (order in thesis)
  dplyr::select(id,
                alpha_if_bg, gamma_if_bg, shape1_if_bg, shape2_if_bg,
                alpha_if_nbg, gamma_if_nbg, shape1_if_nbg, shape2_if_nbg,
                alpha_blue, gamma_blue, shape1_blue, shape2_blue)

save_data(evs.posterior.dep, 
          paste(target_dir, "evs-posterior-dependent-data.rds", sep=FS))

# Generate new dependent tables -------------------------------------------
# get Posterior predictive data
# generate set of dependent tables for each sample from posterior distribution (parameters)
sampled_tables = map_dfr(dep_trials, function(trial_id){
  message(trial_id)
  samples_trial <- posterior_samples.dep %>% filter(id == trial_id) %>% 
    filter(Iteration <= 100) #just use the first 100 samples, not all
  #samples_trial <- evs.posterior.dep %>% filter(id == trial_id) #use evs
  data_trial <- df.dep %>% filter(id == trial_id)
  
  map(seq(1, nrow(samples_trial)), function(idx){
    
    posterior_params <- samples_trial[idx,]
    data_webppl <- list(probs = data_trial,
                        evs_params = posterior_params)

    samples.dep_tables <- webppl(
      program_file = here("webppl-model", "posterior-dependent-data.wppl"),
      data_var = "data",
      model_var = "sample_table",
      data = data_webppl,
      packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS)),
      inference_opts = list(method = "forward", samples = nrow(data_trial))
    ) %>% as_tibble() %>% 
      pivot_wider(names_from = "Parameter", values_from = "value") %>% 
      add_column(id = trial_id)
  })
})
save_data(sampled_tables, 
          paste(target_dir, "sampled-tables-posterior-predictive.rds", sep=FS))

pp_plots_new_tables <- map(dep_trials, function(trial_id){
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
  fn <- paste(target_dir, FS, paste("pp-tables-evs-posterior-", trial_id, 
                                    ".png", sep=""), sep="")
  p <- plot_new_tables(df.trial, sampled_tables, tit, fn)
  return(p)
})

# Posterior predictive ----------------------------------------------------
# Log likelihood plots
ll_fn <- function(df.samples, df.grp){
  message(paste(df.grp$id))
  data_webppl <- list(probs = df.dep %>% filter(id == df.grp$id),
                      samples_posterior = df.samples)
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-dependent-data.wppl"),
    random_seed = params$seed_webppl,
    data_var = "data",
    data = data_webppl,
    packages = c(paste("webppl-model", "node_modules", "dataHelpers", sep = FS))
  ) 
  
  result = tibble(x_blue = samples.pp$x_blue,
                  x_if_bg = samples.pp$x_if_bg,
                  x_if_nbg = samples.pp$x_if_nbg,
                  ll_X_new = samples.pp$ll_X, 
                  ll_obs = samples.pp$ll_obs, 
                  ll_X_blue = samples.pp$ll_x_blue, 
                  ll_X_if_bg = samples.pp$ll_x_if_bg, 
                  ll_X_if_nbg = samples.pp$ll_x_if_nbg, 
                  ll_obs_blue = samples.pp$ll_obs_blue, 
                  ll_obs_if_bg = samples.pp$ll_obs_if_bg, 
                  ll_obs_if_nbg = samples.pp$ll_obs_if_nbg)
  return(result %>% add_column(id = df.grp$id))
}
pp_samples_ll = group_map(posterior_samples.dep %>% group_by(id), ll_fn) %>%
  bind_rows() %>% 
  mutate(trial=id, 
         id = factor(id, levels = c("if1_hh", "if1_uh", "if1_u-Lh", "if1_lh",
                                    "if2_hl", "if2_ul", "if2_u-Ll", "if2_ll"),
                     labels = c(parse(text=expression("if"[1]*":HI")),
                                parse(text=expression(paste(`if`[1], ":UI"))),
                                parse(text=expression("if"[1]*":U"^-{}*"I")),
                                parse(text=expression("if"[1]*":LI")),
                                parse(text=expression("if"[2]*":HL")),
                                parse(text=expression("if"[2]*":UL")),
                                parse(text=expression("if"[2]*":U"^-{}*"L")),
                                parse(text=expression("if"[2]*":LL")))))

# overall log likelihood of data from posterior predictive
ll_X_obs.mean.all = pp_samples_ll %>% group_by(id) %>% 
  summarize(ev = mean(ll_obs), .groups = "keep")
p.all <- pp_samples_ll %>% 
  ggplot(aes(x = ll_X_new)) + geom_density() +
  facet_wrap(~id, scales = "free", ncol = 4, labeller = labeller(id = label_parsed)) +
  geom_point(data = ll_X_obs.mean.all, aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p.all
ggsave(paste(target_dir, "pp-log-likelihood-dependent.png", sep = FS), p.all)

# log likelihood separate for P(b), P(g|b) and P(g|¬b)
pp_samples_ll.long <- pp_samples_ll %>%
  pivot_longer(cols = c(-id, -x_blue, -x_if_bg, -x_if_nbg, -trial), 
               names_to = "prob", names_prefix = "ll_", 
               values_to = "ll")
# ll pp-samples
df.ll_X = pp_samples_ll.long %>% 
  filter(startsWith(prob, "X_") & prob != "X_new") %>% 
  mutate(prob = str_replace(prob, "X_", ""))

# ev ll observed data given params from posterior, seperately for blue/if_bg/if_nbg
ll_X_obs.mean = pp_samples_ll.long %>% group_by(id, trial, prob) %>% 
  filter(startsWith(prob, "obs_")) %>% 
  summarize(ev = mean(ll), .groups = "drop_last") %>% 
  mutate(prob = str_replace(prob, "obs_", ""))

p.if1 <- df.ll_X  %>% filter(startsWith(trial, "if1")) %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(id~prob, ncol=3, scales = "free", 
             labeller = labeller(id = label_parsed, prob = prob_names)) +
  geom_point(data = ll_X_obs.mean %>% filter(startsWith(trial, "if1")), 
             aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p.if1


p.if2 <- df.ll_X %>% filter(startsWith(trial, "if2")) %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(id~prob, ncol=3, scales = "free", 
             labeller = labeller(id = label_parsed, prob = prob_names)) +
  geom_point(data = ll_X_obs.mean %>% filter(startsWith(trial, "if2")), 
             aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p.if2

# returns matrix nb.drawn samples x N=nb.participants, 
# for a single sampled value, given in arg 'col'
format_sampled_ps <- function(df, col) {
  df %>% dplyr::select(idx_sample, id_subj, !!col) %>% 
    group_by(idx_sample) %>% 
    pivot_wider(names_from="id_subj", values_from=col) %>% 
    ungroup() %>% dplyr::select(-idx_sample) %>% as.matrix()
}
# posterior predictive plots (new sampled data for the 3 probabilities)
df.pp_x <- pp_samples_ll %>% dplyr::select(x_blue, x_if_bg, x_if_nbg, id, trial) %>% 
  group_by(id) %>% 
  mutate(idx_sample = seq_along(id)) %>% 
  unnest(cols=c(starts_with("x_"))) %>% 
  mutate(x_if_bg = case_when(is.na(x_if_bg) ~ -1, T ~ x_if_bg),
         x_if_nbg = case_when(is.na(x_if_nbg) ~ -1, T ~ x_if_nbg)) %>% 
  group_by(id, idx_sample) %>% 
  mutate(id_subj = row_number())
         
pp_zoib_plots <- group_map(df.pp_x %>% group_by(trial, id), function(df, grp){
  message(grp$trial)
  df.trial <- df.dep %>% filter(id == grp$trial) %>% 
    mutate(if_bg = case_when(is.na(if_bg) ~ -1, T ~ if_bg),
           if_nbg = case_when(is.na(if_nbg) ~ -1, T ~ if_nbg)) %>% 
    rowid_to_column("idx_subj")
  N = nrow(df.trial)
  y <- df.trial %>% dplyr::select(blue, if_bg, if_nbg) %>%
    as.matrix() %>% matrix(nrow=1, ncol=N*3) %>% as.numeric()
  
  df.mat <- cbind(format_sampled_ps(df, "x_blue"),
                  format_sampled_ps(df, "x_if_bg"),
                  format_sampled_ps(df, "x_if_nbg"))
  p_grps <- c(rep("P(b)", N), rep("P(g|b)", N), rep("P(g|¬b)", N))
  p_grps <- factor(p_grps, levels = c("P(b)", "P(g|b)", "P(g|¬b)"))
  
  
  p <- ppc_dens_overlay_grouped(y = y, yrep = df.mat[1:100,], group=p_grps) + 
    labs(title=trial_names[[grp$trial]])
  
  ggsave(paste(target_dir, paste("pp_zoib_", grp$trial, ".png", sep=""), sep=FS), p)
  return(p)
})


