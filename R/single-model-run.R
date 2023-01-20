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

# use special causal nets, witgh different marginal distributions than uniform
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

# Analyze weights ---------------------------------------------------------

# weights by context
params$likelihoods_zoib <- bind_rows(
  readRDS(here("results/context-free-prior/fit_dep_tbls.rds")), 
  readRDS(here("results/context-free-prior/fit_ind_marginals.rds"))
)
path_par_ind_diffs = here("results/context-free-prior/fit_ind_diffs.rds")
params$likelihoods_gaussian <- readRDS(path_par_ind_diffs) %>% 
  mutate(p = "diff")

weights = run_webppl("webppl-model/weights-contexts.wppl", params) %>% 
  bind_rows(.id = "id") %>% group_by(id) %>% arrange(desc(probs)) %>% 
  unnest(c(support)) %>% 
  unnest(c(table.probs, table.support)) %>% 
  pivot_wider(names_from = table.support, values_from = table.probs) %>% 
  group_by(id) %>% mutate(cdf = cumsum(probs))


# analyze weights
weights %>% group_by(id) %>% summarize(n_states = n())
weights %>% group_by(id) %>% filter(cdf <= 0.99) %>%  summarize(n_states = n())

# compute expected state for each context based on weights 
expected.states = weights %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value") %>% 
  group_by(id, cell) %>%
  summarize(ev_cell = sum(probs * value), .groups = "drop_last") %>% 
  add_column(data = "model states")

# expected slider ratings per context
slider_ratings = data.pe %>% 
  dplyr::select(id, AC, `A-C`, `-AC`, `-A-C`) %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value")
expected.slider_ratings = slider_ratings %>% 
  group_by(id, cell) %>% 
  summarize(ev_cell = mean(value), .groups = "drop_last") %>% 
  add_column(data = "empirical slider ratings")

data.joint <- bind_rows(expected.states, expected.slider_ratings) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))
cell_names <- c(
  `-A-C` = "Neither block falls",
  `-AC` = "Only the green block falls",
  `A-C` = "Only the blue block falls",
  `AC` = "Both blocks fall"
)
# plot expected ratings / model states
plots.expected_states = map(data.joint$id %>% unique, function(id){
  p <- data.joint %>% 
    filter(id == !!id) %>% 
    ggplot(aes(x = cell, y = ev_cell, fill = data, color=data)) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    facet_wrap(~cell, ncol = 2, scales = "free_x", 
               labeller=labeller(cell=cell_names)) +
    theme_minimal() + 
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank()) +
    labs(x = "event", y = "Expected value states/slider ratings", 
         title = id) +
    theme(legend.position = "bottom")
  return(p)
})


#weights by subj and trial
path_model_weights = paste(params$dir_wppl_code, params$fn_get_weights, sep=FS)
weights = run_webppl(path_model_weights, params)

subjs.weights = group_map(data.pe %>% group_by(id), function(df, df.grp){
  trial.id = df.grp$id
  message(trial.id)
  weights.id = weights[[trial.id]][]
  map_dfr(df$prolific_id, function(subj){
    weights.all = weights.id[[subj]] %>% as_tibble() %>% 
      unnest(cols = c(support)) %>% 
      dplyr::select(probs, r, probability, bn_id) %>% 
      add_column(prolific_id = subj, id = trial.id) %>% 
      arrange(desc(probs)) %>% rowid_to_column("idx") %>% 
      mutate(cdf = cumsum(probs))
    # weights = weights.all %>% filter(cdf <= 0.999)
    return(weights.all)
  }) %>% bind_rows()
}) %>% bind_rows() %>% group_by(prolific_id, id) %>% 
  mutate(cn = paste(r, probability, sep = "_")) %>% 
  arrange(desc(probs))

modeled_r = subjs.weights$r %>% unique()
modeled_cn = subjs.weights %>% ungroup() %>% dplyr::select(cn) %>% distinct() %>% pull(cn)
modeled_bn_ids = subjs.weights %>% ungroup() %>% dplyr::select(bn_id) %>%
  distinct() %>% pull(bn_id)

data_pe.weights = left_join(
  data.pe %>% ungroup() %>% 
    dplyr::select(prolific_id, id, AC, `A-C`, `-AC`, `-A-C`) %>% 
    # mutate(cn = list(modeled_cn)) %>% unnest(cols = c(cn)),
    mutate(bn_id = list(modeled_bn_ids)) %>% unnest(cols = c(bn_id)),
  subjs.weights %>% dplyr::select(prolific_id, id, bn_id, probs)
) %>% 
  mutate(probs = case_when(is.na(probs) ~ 0, T ~ probs)) %>% 
  mutate(bn_id.tmp = bn_id) %>% 
  separate(bn_id.tmp, into = c("r", "probability", "t1", "t2", "t3", "t4"), sep = "_") %>% 
  dplyr::select(-t1, -t2, -t3, -t4)

data.weights = left_join(
  data_pe.weights, 
  data.uc %>% dplyr::select(prolific_id, id, utterance)
) %>% rename(uc_task = utterance)

# cn: r + probability
p_cn_Dij = data.weights %>% 
  group_by(prolific_id, id, r, AC, `A-C`, `-AC`, `-A-C`, uc_task, probability) %>% 
  summarize(p = sum(probs), .groups = "drop_last") %>% 
  mutate(r = case_when(r == "A implies C" ~ "A -> C", 
                       r == "-A implies C" ~ "-A -> C",
                       r == "C implies A" ~ "C -> A",
                       r == "-C implies A" ~ "- C-> A",
                       r == "A || C" ~ "ind")) %>% 
  arrange(desc(p))
p_r_Dij = p_cn_Dij %>% summarize(p = sum(p), .groups = "drop") %>% 
  group_by(prolific_id, id) %>% 
  arrange(desc(p)) 

# plot posterior P(r|D_ij) for all subjects, mark average posterior P(r|D_ij) 
# for each context i
p_cn_Di.avg = p_cn_Dij %>% group_by(id, r, probability) %>% 
  summarize(p = mean(p), .groups = "drop_last") %>% 
  arrange(desc(p))
p_r_Di.avg = p_r_Dij %>% group_by(id, r) %>% 
  summarize(p = mean(p), .groups = "drop_last")

plot.posterior_r = p_r_Dij %>%
  ggplot(aes(x = r, y = p)) +
  geom_point(aes(group = prolific_id)) + geom_line(aes(group = prolific_id)) + 
  facet_wrap(~id, scales = "free_x") +
  geom_point(data = p_r_Di.avg, color = "red")
plot.posterior_r

target_dir = here(params$dir_results)
ggsave(paste(target_dir, "p_r_Dij.png", sep=FS), plot.posterior_r, 
       width = 10, height = 6)
# Run Model ---------------------------------------------------------------
path_model_file = paste(params$dir_wppl_code, params$fn_rsa_single_run, sep=FS)
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

plot_correlation(freq.joint, color = "id") + facet_wrap(~id)
plot_correlation(freq.joint.types, color = "id") + facet_wrap(~id)

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








