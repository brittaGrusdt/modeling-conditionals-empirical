library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggpubr)

source(here("R", "generate-model-states.R"))
source(here("R", "plot-functions.R"))
source(here("R", "prediction-functions.R"))
# Runs default-model
# Setup -------------------------------------------------------------------
# Prior
# active_config = "cf_prior_match_tables"
active_config = "cf_prior_match_kl"
# active_config = "context_sensitive_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)

# Data
Sys.setenv(R_CONFIG_ACTIVE = "cleaned_data")
params.data = config::get()
path_empiric_tbls_ids = paste(params.data$dir_data, 
                              params.data$fn_tbls_empiric_pids, sep = FS)

# Behavioral data
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task.smooth, slider) %>% 
  translate_standardized2model() 
data.uc = data.behav %>%  filter(!is.na(uc_task)) %>% 
  dplyr::select(prolific_id, id, utterance) %>% 
  add_column(n = 1) %>% group_by(id, utterance)
data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task.smooth) %>% 
  pivot_wider(names_from = "utt.standardized", values_from = "pe_task.smooth") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`)

data.observed = left_join(data.uc, data.pe)
empirical.context = data.observed %>% group_by(id, utterance) %>% 
  summarize(n = n(), .groups = "drop_last") %>% 
  mutate(behavioral = n / sum(n))

# set states if presampled by running prior once before running model
dir_presampled_states = here(params$dir_data, params$name, 
                             params$dir_presampled_states)
if(!dir.exists(dir_presampled_states)) {
  dir.create(dir_presampled_states, recursive = T)
} 
fn_presampled_states = paste("states-", params$seed_webppl, ".rds", sep="")
fn_kl_divs =  paste("KL-divergences-", params$seed_webppl, ".rds", sep="")
path_presampled = paste(dir_presampled_states, fn_presampled_states, sep = FS)
path_kl_divs <- paste(dir_presampled_states, fn_kl_divs, sep = FS)

if(params$states_presampled){
  if(!file.exists(path_presampled)) {
    states = generate_model_states(active_config)
    save_data(states, path_presampled)
  } else {
    states <- readRDS(path_presampled)
  }
  params$prior_samples <- states
}

# Run Model ---------------------------------------------------------------
path_model_file = here(params$dir_webppl_model, params$fn_rsa_model)
params$packages <- params$packages[1:2]
params$p_utts = rep(1 / length(params$utterances), length(params$utterances))
posterior <- run_webppl(path_model_file, params)
# restructure data and save
speaker <- posterior$distributions %>% 
  structure_speaker_data(params) %>% group_by(bn_id)

if(params$states_presampled) {
  tables.model = states %>% unnest(c(table.probs, table.support)) %>% 
    group_by(bn_id) %>% 
    pivot_wider(names_from = "table.support", values_from = "table.probs")
} else {
  tables.model <- posterior$distributions$bns %>% as_tibble() %>% 
    dplyr::select(bn_id, table.probs, table.support) %>% 
    unnest(c(table.probs, table.support)) %>% 
    group_by(bn_id) %>% 
    pivot_wider(names_from = "table.support", values_from = "table.probs")
}

# Model Predictions -------------------------------------------------------
# --- Predictions by Relation (for Dirichlet) ---
if(active_config == "context_sensitive_prior") {
  # for context-sensitive prior, matching is not needed: the tables were sampled
  # from Dirichlet distributions fitted to participants' data for each context
  predictions.context = predictions_by_relation(speaker) %>% 
    rename(id = r)

} else if(active_config == "cf_prior_match_tables") {
  # --- Predictions by matching tables ---
  df.predictions = predictions_by_matching_tables(speaker, tables.model, 
                                                  path_empiric_tbls_ids)
  predictions.context = df.predictions$predictions
  df.matches <- df.predictions$matches
  # Check matched tables: how many trials were matched? (i.e.taken into account)
  n_trials.match = df.matches %>% group_by(subject_id) %>% dplyr::count() %>% 
    arrange(desc(n)) %>% pull(n) %>% sum()
  n_trials.all = nrow(data.uc)
  print(n_trials.match / n_trials.all)

} else if(active_config == "cf_prior_match_kl") {
  # --- Predictions based on states with lowest KL-divergences ---
  tbls.empiric <- data.pe %>% 
    dplyr::select(prolific_id, id, AC, `A-C`, `-AC`, `-A-C`)
  if(!file.exists(path_kl_divs)) {
    kl.divergences = compute_kl_divergences(tables.model, tbls.empiric)
    save_data(kl.divergences, path_kl_divs)
  } else {
    kl.divergences <- readRDS(path_kl_divs)
  } 
  
  predictions.context = left_join(kl.divergences %>% filter(idx == 1), 
                                  speaker) %>% 
    group_by(id, utterance) %>% 
    summarize(model = mean(probs), .groups = "drop_last") %>% 
    arrange(desc(model)) %>% 
    mutate(utt = utterance) %>% 
    chunk_utterances() %>% 
    rename(utt_type = utterance, utterance = utt)
  
  
  tbls.most_similar = left_join(
    data.pe %>% dplyr::select(c(prolific_id, id, AC, `A-C`, `-AC`, `-A-C`)),
    kl.divergences %>% 
      dplyr::select(c(prolific_id, id, AC.match, `A-C.match`, 
                      `-AC.match`, `-A-C.match`, `kl.div`, idx))
  ) %>% 
    mutate(AC = round(AC, 2), 
           `A-C` = round(`A-C`, 2), 
           `-AC` = round(`-AC`, 2), 
           `-A-C` = round(`-A-C`, 2),
           AC.match = round(as.numeric(AC.match), 2), 
           `A-C.match` = round(as.numeric(`A-C.match`), 2), 
           `-AC.match` = round(as.numeric(`-AC.match`), 2), 
           `-A-C.match` = round(as.numeric(`-A-C.match`), 2))
  
  # worst matches of best matches:
  tbls.most_similar %>% filter(idx == 1) %>% arrange(desc(kl.div))
  tbls.most_similar %>% filter(idx == 1) %>% pull(kl.div) %>% summary()
  
  p <- tbls.most_similar %>% 
    mutate(states = case_when(idx == 1 ~ "most similar",
                              T ~ "other")) %>% 
    ggplot(aes(x=kl.div)) + geom_density(aes(color = states)) +
    labs(x = "KL-divergence") + theme(legend.position = "top")
}



# Joint data empirical + model --------------------------------------------
results.joint <- left_join(predictions.context, empirical.context) %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)),
         stimulus = id) %>% 
  separate("stimulus", into = c("relation", "prior"), sep = "_") %>% 
  dplyr::select(-prior)
  
results.joint.utt_type <- results.joint %>% 
  group_by(id, relation, utt_type) %>% 
  summarize(model = sum(model), behavioral = sum(behavioral),
            .groups = "keep")


# Results -----------------------------------------------------------------
predictions.context %>% filter(id == "if1_uh")
predictions.context %>% filter(id == "independent_uh")

# plots
grp_var = "utterance"
p_scatter.types = plot_correlation(
  results.joint.utt_type %>% rename(utterance = utt_type),
  color = grp_var
)
p_scatter.types

p_scatter.utts = plot_correlation(results.joint, label.y = 0.3, color=grp_var) 
p_scatter.utts + facet_wrap(as.formula(paste("~", grp_var)))

