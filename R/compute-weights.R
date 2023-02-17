library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "plot-functions.R"))

cell_names <- c(
  `-A-C` = "Neither block falls",
  `-AC` = "Only the green block falls",
  `A-C` = "Only the blue block falls",
  `AC` = "Both blocks fall"
)
colors <- c(`empirical slider ratings` = "darkgreen", `model states` = "firebrick")
# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
result_dir = params$dir_results


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

slider_ratings = data.pe %>% 
  dplyr::select(id, prolific_id, AC, `A-C`, `-AC`, `-A-C`) %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value")

data.observed = left_join(data.uc, data.pe)
params$observations = data.observed

# draw samples from prior
prior <- webppl(program_file = here("webppl-model", "run-state-prior.wppl"),
                data_var = "data",
                data = params,
                random_seed = params$seed_webppl,
                packages = params$packages[1:2]) 
states = prior$prior$support %>% as_tibble()
params$prior_samples = states
# params$p_utts = rep(1 / length(params$utterances), length(params$utterances))

# load likelihood parameters (for fitted data)
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

# Weights by contexts (M_Ci) ----------------------------------------------
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))
weights_ci = run_webppl("webppl-model/weights-contexts.wppl", params) %>% 
  bind_rows(.id = "id") %>% group_by(id) %>% arrange(desc(probs)) %>% 
  unnest(c(support)) %>% 
  unnest(c(table.probs, table.support)) %>% 
  pivot_wider(names_from = table.support, values_from = table.probs) %>% 
  group_by(id) %>% mutate(cdf = cumsum(probs))


save_data(weights_ci %>% dplyr::select(id, bn_id, probs), 
          paste(result_dir, "weights_ci.rds", sep=FS))


# analyze weights
weights_ci %>% group_by(id) %>% summarize(n_states = n())
weights_ci %>% group_by(id) %>% filter(cdf <= 0.99) %>%  summarize(n_states = n())

# relations
weights_ci.relations = weights_ci %>% group_by(id, r) %>%
  summarize(p = sum(probs), .groups = "drop_last") %>% 
  arrange(id, desc(p))
weights_ci.relations %>% filter(id == "independent_ul")
weights_ci.relations %>% filter(startsWith(id, "if1"))
weights_ci.relations %>% filter(startsWith(id, "if2"))


# compute expected state for each context based on weights 
expected.states.ci = weights_ci %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value") %>% 
  group_by(id, cell) %>%
  summarize(ev_cell = sum(probs * value), .groups = "drop_last") %>% 
  add_column(data = "model states")

# expected slider ratings per context
expected.slider_ratings.ci = slider_ratings %>% 
  group_by(id, cell) %>% 
  summarize(ev_cell = mean(value), .groups = "drop_last") %>% 
  add_column(data = "empirical slider ratings")

data.joint.ci <- bind_rows(expected.states.ci, expected.slider_ratings.ci) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))

# plot expected ratings / model states
plots.expected_states.ci = map(data.joint.ci$id %>% unique, function(id){
  
  p <- data.joint.ci %>% 
    filter(id == !!id) %>% 
    ggplot(aes(x = cell, fill = data, color=data)) +
    geom_bar(aes(y = ev_cell), stat = "identity", 
             position = position_dodge()) + 
    facet_wrap(~cell, ncol = 2, scales = "free_x", 
               labeller=labeller(cell=cell_names)) +
    theme_minimal() + 
    scale_color_manual(name = "data", values = colors) +
    scale_fill_manual(name = "data", values = colors) +
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank()) +
    labs(x = "event", y = "Expected value states/slider ratings", 
         title = id) +
    theme(legend.position = "bottom")
  return(p)
})
plots.expected_states.ci


# Weights by subj and trial (M_dij) ---------------------------------------
weights_dij = run_webppl("webppl-model/weights-dij.wppl", params) 
subjs.weights = group_map(data.pe %>% group_by(id), function(df, df.grp){
  trial.id = df.grp$id
  message(trial.id)
  weights_dij.id = weights_dij[[trial.id]][]
  map_dfr(df$prolific_id, function(subj){
    weights_dij.all = weights_dij.id[[subj]] %>% as_tibble() %>% 
      unnest(cols = c(support)) %>% 
      dplyr::select(probs, r, probability, bn_id) %>% 
      add_column(prolific_id = subj, id = trial.id) %>% 
      arrange(desc(probs)) %>% rowid_to_column("idx") %>% 
      mutate(cdf = cumsum(probs))
    return(weights_dij.all)
  }) %>% bind_rows()
}) %>% bind_rows() %>% group_by(prolific_id, id) %>% 
  mutate(cn = paste(r, probability, sep = "_")) %>% 
  arrange(desc(probs))

save_data(subjs.weights %>% dplyr::select(prolific_id, id, bn_id, probs), 
          paste(result_dir, "weights_dij.rds", sep=FS))

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

data_pe_uc.weights.dij = left_join(
  data_pe.weights, 
  data.uc %>% dplyr::select(prolific_id, id, utterance)
) %>% rename(uc_task = utterance) %>% 
  arrange(prolific_id, id, desc(probs)) %>% 
  group_by(prolific_id, id) %>% mutate(cdf = cumsum(probs))

# look at weights for specific context and subject
# and compute expected state for each context and subject based on weights 
expected.slider_ratings.dij = slider_ratings %>% 
  add_column(data = "empirical slider ratings") %>% 
  rename(ev_cell = value)
expected.states.dij = data_pe_uc.weights.dij %>% 
  dplyr::select(prolific_id, id, bn_id, probs) %>% 
  separate("bn_id", into = c("r", "probabilities", "AC", "A-C", "-AC", "-A-C"), sep = "_") %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value") %>% 
  mutate(value = as.numeric(value)) %>% 
  group_by(prolific_id, id, cell) %>%
  summarize(ev_cell = sum(probs * value), .groups = "drop_last") %>% 
  add_column(data = "model states")

expected.joint.dij <- bind_rows(expected.states.dij, expected.slider_ratings.dij) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))

context <- "if1_uh"
subj1 <- "5f68d4c59b202b055cc873e5"
subj2 <- "5f04c46bcfd81678efd229de"
expected.joint.dij %>%
  filter(id == !!context & (prolific_id == !!subj1)) %>% 
  ggplot(aes(x = cell, fill = data, color=data)) +
  geom_bar(aes(y = ev_cell), stat = "identity", 
           position = position_dodge()) + 
  facet_wrap(~cell, ncol = 2, scales = "free_x", 
             labeller=labeller(cell=cell_names)) +
  scale_color_manual(name = "data", values = colors) +
  scale_fill_manual(name = "data", values = colors) +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  labs(x = "event", y = "Expected value states/slider ratings")



# weights by most similar states M_sim -------------------------------------
weights_sim = run_webppl("webppl-model/weights-sim.wppl", params) 
subjs.weights.sim = group_map(data.pe %>% group_by(id), function(df, df.grp){
  trial.id = df.grp$id
  message(trial.id)
  weights_sim.id = weights_sim[[trial.id]][]
  map_dfr(df$prolific_id, function(subj){
    weights_sim.all = weights_sim.id[[subj]] %>% as_tibble() %>% 
      unnest(cols = c(support)) %>% 
      dplyr::select(probs, r, probability, bn_id) %>% 
      add_column(prolific_id = subj, id = trial.id) %>% 
      arrange(desc(probs)) %>% rowid_to_column("idx") %>% 
      mutate(cdf = cumsum(probs))
    return(weights_sim.all)
  }) %>% bind_rows()
}) %>% bind_rows() %>% group_by(prolific_id, id) %>% 
  mutate(cn = paste(r, probability, sep = "_")) %>% 
  arrange(desc(probs))

save_data(subjs.weights.sim %>% dplyr::select(prolific_id, id, bn_id, probs), 
          paste(result_dir, "weights_dij_sim.rds", sep=FS))


expected.states.dij.sim = subjs.weights.sim %>% 
  dplyr::select(prolific_id, id, bn_id) %>% 
  separate("bn_id", into = c("r", "probabilities", "AC", "A-C", "-AC", "-A-C"), sep = "_") %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value") %>% 
  mutate(value = as.numeric(value)) %>% 
  add_column(data = "model states") %>% 
  dplyr::select(-r, -probabilities)

expected.joint.dij.sim <- bind_rows(
  expected.states.dij.sim, 
  expected.slider_ratings.dij %>% rename(value = ev_cell)
) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))


context <- "if1_uh"
subj1 <- "5f68d4c59b202b055cc873e5"
subj2 <- "5f04c46bcfd81678efd229de"
expected.joint.dij.sim %>%
  filter(id == !!context & (prolific_id == !!subj2)) %>% 
  ggplot(aes(x = cell, fill = data, color=data)) +
  geom_bar(aes(y = value), stat = "identity", 
           position = position_dodge()) + 
  facet_wrap(~cell, ncol = 2, scales = "free_x", 
             labeller=labeller(cell=cell_names)) +
  scale_color_manual(name = "data", values = colors) +
  scale_fill_manual(name = "data", values = colors) +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  labs(x = "event", y = "Expected value states/slider ratings")

# Plots  ------------------------------------------------------------------
# Plot distribution over relations for each participant for weights computed
# by participant and trial (M_dij)
p_cn_Dij = data_pe_uc.weights.dij %>% 
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

