library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggpubr)

# Runs default-model
# Setup -------------------------------------------------------------------
# Priors
# active_config = "situation_specific_prior"
active_config = "abstract_default_params"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)

# Data
Sys.setenv(R_CONFIG_ACTIVE = "cleaned_data")
params.data = config::get()
path_empiric_tbls_ids = paste(params.data$dir_data, 
                              params.data$fn_tbls_empiric_pids, sep = FS)
                              
cns = create_causal_nets(params$p_noise, params$p_cp, params$p_ant, 
                         params$rels_dep, params$p_a, params$p_c)
params$causal_nets_dep = cns$dep
params$causal_nets_ind = cns$ind

## Generate/Retrieve utterances
path_utterances = here(params$dir_model_input, params$fn_utterances)
if(params$generate_utterances || !file.exists(path_utterances)){
  utterances <- generate_utts(params, path_to_target = path_utterances)
} else {
  utterances <- readRDS(path_utterances)
  print(paste("utterances read from:", path_utterances))
}
params$utterances <- utterances


prior <- webppl(program_file = here("webppl-model", "run-state-prior.wppl"),
                data_var = "data",
                data = params,
                packages = params$packages[1:2]) 
states = prior$prior$support %>% as_tibble()
params$prior_samples = states

tables = states %>% unnest(c(table.probs, table.support)) %>% 
  group_by(bn_id) %>% 
  pivot_wider(names_from = "table.support", values_from = "table.probs") %>% 
  mutate(AC.round=as.integer(round(AC, 2) * 100),
         `A-C.round`=as.integer(round(`A-C`,2) * 100),
         `-AC.round`=as.integer(round(`-AC`, 2) * 100),
         `-A-C.round`=as.integer(round(`-A-C`, 2) * 100)) %>% 
  rowid_to_column("table_id")  
tables.empiric_pids = match_sampled_and_empiric_tables(tables, 
                                                       path_empiric_tbls_ids) %>% 
  ungroup() %>% 
  dplyr::select(r, probability, bn_id, empirical_id, stimulus, p_id, match.empirical)

# Run Model ---------------------------------------------------------------
path_model_file = here(params$dir_webppl_model, params$fn_rsa_model)
params$packages <- params$packages[1:2]
params$p_utts = rep(1 / length(params$utterances), length(params$utterances))
posterior <- run_webppl(path_model_file, params)

# restructure data and save
speaker <- posterior$distributions %>% 
  structure_speaker_data(params) %>% group_by(bn_id)


# Predictions by matching tables ------------------------------------------
predictions <- left_join(speaker, tables.empiric_pids) %>%
  filter(match.empirical) %>% unnest(c(p_id)) %>% 
  dplyr::select(bn_id, utterance, probs, p_id) %>% 
  separate("p_id", into = c("subject_id", "r", "prior"), sep = "_") %>% 
  unite("id", "r", "prior", sep = "_")

df.matches = predictions %>% dplyr::select(bn_id, subject_id, id) %>% 
  distinct() %>% group_by(subject_id, id) %>% dplyr::count() %>% arrange(desc(n))
df.matches %>% filter(n>4)

predictions.context = predictions %>% group_by(utterance, id) %>% 
  summarize(model = mean(probs), .groups = "drop_last") %>% 
  arrange(desc(model)) %>% 
  mutate(utt = utterance) %>% 
  chunk_utterances() %>% rename(utt_type = utterance, utterance = utt)

predictions.context %>% filter(id == "if1_uh")
predictions.context %>% filter(id == "independent_uh")


# Behavioral data ---------------------------------------------------------
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
  mutate(behavioral.p = n / sum(n))

results.joint <- left_join(predictions.context, empirical.context) %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
results.joint.utt_type <- results.joint %>% group_by(id, utt_type) %>% 
  summarize(model = mean(model), behavioral.p = mean(behavioral.p))

p_scatter.types = results.joint.utt_type %>% 
  ggscatter(x = "behavioral.p", y = "model", add = "reg.line",
            conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
            xlab = "Empirical observations", ylab = "Model predictions") +
  geom_point(size=1.5, aes_string(x="behavioral.p", y="model", 
                                  color="utt_type"))
p_scatter.types
p_scatter.types + facet_wrap(~id)

p_scatter.utts = results.joint %>% 
  ggscatter(x = "behavioral.p", y = "model", add = "reg.line",
            conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
            xlab = "Empirical observations", ylab = "Model predictions") +
  geom_point(size=1.5, aes_string(x="behavioral.p", y="model", color="utterance"))
p_scatter.utts
p_scatter.utts + facet_wrap(~id)

