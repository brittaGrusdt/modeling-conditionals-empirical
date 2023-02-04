library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "plot-functions.R"))

theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)

dep_data_dir = paste(here(params$dir_results), "dependent-contexts", sep=FS)
ind_data_dir = paste(here(params$dir_results), "independent-contexts", sep=FS)

# Data
Sys.setenv(R_CONFIG_ACTIVE = "cleaned_data")
params.data = config::get()

# use special causal nets, with different marginal distributions than uniform
# cns = create_causal_nets(params$p_noise, params$p_cp, params$p_ant,
#                          params$rels_dep, params$p_a, params$p_c)
# params$causal_nets_dep = cns$dep
# params$causal_nets_ind = cns$ind

# behavioral data
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
  translate_standardized2model() 

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

data.observed = left_join(data.uc, data.pe)
params$observations = data.observed


# draw samples from prior
prior <- webppl(program_file = here("webppl-model", "run-state-prior.wppl"),
                data_var = "data",
                data = params,
                packages = params$packages[1:2]) 
states = prior$prior$support %>% as_tibble()
params$prior_samples = states
params$p_utts = rep(1 / length(params$utterances), length(params$utterances))

# load likelihood parameters fitted to data
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

params$likelihoods_zoib <- bind_rows(pars.dep, pars.ind)
params$likelihoods_bernoulli <- pars.ind.theta_p

# Run Model ---------------------------------------------------------------
# 1. predictions by contexts
path_model_file = paste(params$dir_wppl_code, params$fn_rsa_single_run, sep=FS)

# include observations as ratios and counts
trials <- data.pe$id %>% unique()
params$observed_utts_ratios <- left_join(
    tibble(id = trials, utterance = list(params$utterances)) %>%
      unnest(c(utterance)),
    params$observations %>% group_by(id, utterance) %>% count(), 
  ) %>% mutate_if(is.numeric, coalesce, 0) %>% 
  group_by(id) %>% 
  mutate(ratio = n / sum(n))

params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))
posterior <- run_webppl(path_model_file, params)
model_predictions.ci <- posterior %>% 
  map(function(x){as_tibble(x) %>% mutate(ll_ci = as.numeric(ll_ci))}) %>% 
  bind_rows() %>% group_by(id) %>% 
  mutate(p_hat_round = round(p_hat, 2)) %>% 
  arrange(desc(p_hat))

production.joint.ci = left_join(model_predictions.ci, 
                             params$observed_utts_ratios %>% rename(p = ratio)) %>% 
  rename(model = p_hat, behavioral = p) %>% 
  mutate(relation = case_when(startsWith(id, "independent") ~ "independent", 
                              startsWith(id, "if1") ~ "if1",
                              startsWith(id, "if2") ~ "if2"))

plot_correlation(production.joint.ci, color = "id") + facet_wrap(~id)




# Predictions by each participant -----------------------------------------
path_model_file = paste(params$dir_wppl_code, "model-single-run-by-dij.wppl", sep=FS)
posterior <- run_webppl(path_model_file, params)
model_predictions.dij <- posterior %>% 
  map(function(x){
    as_tibble(x$predictions) %>% 
      mutate(ll_subj = as.numeric(ll_subj), 
             id = x$id)
    }) %>% 
  bind_rows() %>% group_by(id, prolific_id) %>% 
  unnest(c(prediction.probs, prediction.support)) %>% 
  rename(p_hat = prediction.probs, utterance = prediction.support) %>% 
  mutate(p_hat_round = round(p_hat, 2)) %>% 
  arrange(desc(p_hat))

# prediction (1xnb_utterances) for each participant and trial
production.joint.dij = left_join(
  model_predictions.dij, 
  params$observations %>% dplyr::select(prolific_id, id, utterance) %>% 
    rename(uc_task = utterance)
) %>% rename(model = p_hat)

# average prediction across all participants
model.avg <- production.joint.dij %>% group_by(id, utterance) %>% 
  summarize(model = mean(model), .groups = "drop_last")
behav.avg <- data.uc %>% 
  group_by(id, utterance) %>%
  summarize(behavioral = n(), .groups = "drop_last") %>% 
  mutate(behavioral = behavioral / sum(behavioral))
joint.avg.dij <- left_join(model.avg, behav.avg) %>% 
  mutate(behavioral = case_when(is.na(behavioral) ~ 0, T ~ behavioral))

joint.avg.dij %>% 
  ggplot(aes(x=behavioral, y = model)) +
  geom_point(aes(color = utterance)) +
  geom_smooth(method = 'lm')

joint.avg.dij %>% 
  ggplot(aes(x=behavioral, y = model)) +
  geom_point(aes(color = utterance)) +
  facet_wrap(~id)

# just models best prediction per dij
model.best <- production.joint.dij %>% filter(model == max(model)) %>% 
  group_by(id, utterance) %>% 
  summarize(n = n(), .groups = "drop_last") %>% 
  mutate(model = n / sum(n)) %>% dplyr::select(-n)

# plot for single participant
production.joint.dij %>% filter(prolific_id == "5ee454ad40c72e152991826e") %>% 
  mutate(utterance = factor(utterance, 
                            levels = c("might A", "might C", "might -A", "might -C",
                                      "C and A", "C and -A", "-C and A", "-C and -A",
                                      "A", "C", "-A", "-C",
                                      "A > C", "A > -C", "-A > C", "-A > -C",
                                      "C > A", "C > -A", "-C > A", "-C > -A"
  ))) %>% 
  ggplot(aes(x = model, y = utterance)) +
  geom_bar(stat = "identity") + 
  geom_point(aes(x = 0, y = uc_task), size=1.5, color = 'red') + 
  labs(y = "utterance", x = "predicted speaker probability") + 
  facet_wrap(~id, scales = "free_y")









