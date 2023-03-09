library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggthemes)
source(here("R", "herlpers-plotting.R"))
source(here("R", "helpers-load-data.R"))

theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))

cell_names <- c(
  `-A-C` = "¬b¬g",
  `-AC` = "¬bg",
  `A-C` = "b¬g",
  `AC` = "bg"
)
colors <- c(`empirical slider ratings` = "darkgreen", `model states` = "firebrick")
# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
result_dir = params$dir_results

dep_data_dir = paste(here(params$dir_results), "dependent-contexts", "zoib-model", sep=FS)
ind_data_dir = paste(here(params$dir_results), "independent-contexts", "zoib-model", sep=FS)

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
pars_likelihoods <- get_likelihood_params_fitted_data(params)
params$likelihoods_zoib <- pars_likelihoods$likelihoods_zoib
params$likelihoods_gaussian <- pars_likelihoods$likelihoods_gaussian

# Weights by contexts (M_Ci) ----------------------------------------------
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()

weights = map(c("uninformative", "informative"), function(prior_r){
  params$prior_relations <- par_relations[[prior_r]]
  weights_ci = run_webppl("webppl-model/weights-contexts.wppl", params) %>% 
    bind_rows(.id = "id") %>% group_by(id) %>% arrange(desc(probs)) %>% 
    unnest(c(support)) %>% 
    unnest(c(table.probs, table.support)) %>% 
    pivot_wider(names_from = table.support, values_from = table.probs) %>% 
    group_by(id) %>% mutate(cdf = cumsum(probs)) %>% 
    add_column(prior_r = prior_r)
}) %>% bind_rows() 

save_data(weights %>% dplyr::select(id, bn_id, probs, prior_r), 
          paste(result_dir, FS, "weights_ci_", params$n_forward_samples, ".rds", 
                sep=""))

# analyze weights
weights %>% group_by(prior_r, id) %>% summarize(n_states = n())
weights %>% group_by(prior_r, id) %>% filter(cdf <= 0.99) %>%
  summarize(n_states = n())

# relations
weights.relations = weights %>% filter(prior_r == "uninformative") %>% 
  group_by(id, r) %>% summarize(p = sum(probs), .groups = "drop_last") %>% 
  arrange(id, desc(p)) %>% mutate(p = round(p, 2)) %>% 
  filter(p != 0)
weights.relations %>% filter(id == "independent_ul")
weights.relations %>% filter(startsWith(id, "if1"))
weights.relations %>% filter(startsWith(id, "if2"))


# compute expected state for each context based on weights 
expected.states = weights %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value") %>% 
  group_by(prior_r, id, cell) %>%
  summarize(ev_cell = sum(probs * value), .groups = "drop_last") %>% 
  add_column(data = "model states")

# expected slider ratings per context
expected.slider_ratings = slider_ratings %>% 
  group_by(id, cell) %>% 
  summarize(ev_cell = mean(value), .groups = "drop_last") %>% 
  add_column(data = "empirical slider ratings")

data.joint <- bind_rows(
  expected.states, 
  expected.slider_ratings %>% add_column(prior_r = "informative"),
  expected.slider_ratings %>% add_column(prior_r = "uninformative")
) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))

# plot expected ratings / model states
plots.expected_states = map(data.joint$id %>% unique, function(id){
  p <- data.joint %>% 
    filter(id == !!id) %>% 
    ggplot(aes(x = cell, fill = data, color=data)) +
    geom_bar(aes(y = ev_cell), stat = "identity", 
             position = position_dodge()) + 
    facet_grid(prior_r~cell, scales = "free_x", 
               labeller=labeller(cell=cell_names)) +
    scale_color_manual(name = "data", values = colors) +
    scale_fill_manual(name = "data", values = colors) +
    theme(axis.ticks.x = element_blank()) + 
    labs(x = "event", y = "Expected value states/slider ratings", title = id) +
    theme(legend.position = "bottom")
  return(p)
})
plots.expected_states




