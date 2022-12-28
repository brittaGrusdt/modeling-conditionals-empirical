library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "plot-functions.R"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "context_free_prior"
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

# single model run
path_model_file = paste(params$dir_wppl_code, params$fn_rsa_single_run, sep=FS)
params$p_utts = rep(1 / length(params$utterances), length(params$utterances))

posterior <- run_webppl(path_model_file, params)

model_predictions <- posterior %>% 
  map(function(x){as_tibble(x)}) %>% bind_rows() %>% 
  unnest(cols = c(predictions)) %>% 
  rename(p_hat = prediction.probs, utterance = prediction.support) 

data.joint = left_join(model_predictions, 
                       data.uc %>% dplyr::select(-n) %>% 
                         rename(observed = utterance), 
                       by=c("id", "prolific_id")) %>% 
  unnest(cols = c(p_hat, utterance)) %>% group_by(prolific_id, id) %>% 
  mutate(predicted = utterance) %>% chunk_utterances() %>% 
  rename(predicted.type = utterance) %>% 
  mutate(utterance = observed) %>% chunk_utterances() %>% 
  rename(observed.type = utterance)

frequencies.predicted = data.joint %>% group_by(id, predicted) %>% 
  summarize(model = mean(p_hat), .groups = "drop_last") %>% 
  rename(utterance = predicted) %>% arrange(desc(model))
frequencies.observed = data.joint %>% 
  dplyr::select(prolific_id, id, observed, observed.type) %>% distinct() %>% 
  group_by(id, observed) %>% 
  summarize(n = n(), .groups = "drop_last") %>% 
  mutate(behavioral = n / sum(n)) %>% rename(utterance = observed)

freq.joint <- left_join(frequencies.predicted, frequencies.observed, 
                        by = c("id", "utterance")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)),
         utt = utterance, stimulus = id) %>% 
  chunk_utterances() %>% rename(utt_type = utterance, utterance = utt) %>% 
  separate("stimulus", into = c("relation", "prior"), sep = "_") %>% 
  dplyr::select(-prior)

# freq.joint.long <- freq.joint %>% dplyr::select(-n) %>% 
#   pivot_longer(cols = c(model, behavioral), names_to="predictor", values_to="p") 

freq.joint.types = freq.joint %>% group_by(id, relation, utt_type) %>% 
  summarize(model = sum(model), behavioral = sum(behavioral),
            .groups = "drop_last") %>% rename(utterance = utt_type)
# Some plots --------------------------------------------------------------
p = plot_correlation(freq.joint)
p + facet_wrap(~utterance)

p.types = plot_correlation(freq.joint.types)
p.types




# most likely prediction = observed utterance / utterance type?
df.most_likely_prediction = data.joint %>% filter(p_hat == max(p_hat)) %>% 
  mutate(correct_utt = predicted == observed, 
         correct_type = predicted.type == observed.type)

df.most_likely_prediction %>% group_by(id) %>% 
  summarize(correct_utt = mean(correct_utt), 
            correct_type = mean(correct_type)) %>% 
  arrange(desc(correct_type))

# plot for single participant
data.joint %>% filter(prolific_id == "5ee454ad40c72e152991826e") %>% 
  mutate(predicted = factor(predicted, levels = c("might A", "might C", "might -A", "might -C",
                                                  "C and A", "C and -A", "-C and A", "-C and -A",
                                                  "A", "C", "-A", "-C",
                                                  "A > C", "A > -C", "-A > C", "-A > -C",
                                                  "C > A", "C > -A", "-C > A", "-C > -A"
  ))) %>% 
  # ggplot(aes(y = reorder_within(utt, p_hat, id), x = p_hat)) + 
  ggplot(aes(x = p_hat, y = predicted)) +
  geom_bar(stat = "identity") + 
  geom_point(aes(x = 0, y = observed), size=1.5, color = 'red') + 
  #scale_y_reordered() +
  labs(y = "utterance", x = "predicted speaker probability") + 
  facet_wrap(~id, scales = "free_y")








