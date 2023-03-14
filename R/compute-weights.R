library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggthemes)
source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))

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

################################################################################
# use more fine-grained causal nets (causal power, noise: low/unc/high)
active_config = "fine_grained_causal_nets"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
pars.cns <- config::get()
cns = create_causal_nets(pars.cns$dep_causal_power, pars.cns$dep_noise, 
                         pars.cns$dep_marginal, params$rels_dep,
                         pars.cns$ind_a, pars.cns$ind_c)
# dont allow causal power to be lower category than noise
params$causal_nets_dep = cns$dep[!str_detect(cns$dep, "-low-") &
                                   !str_detect(cns$dep, "-unc-high")]
params$causal_nets_ind = cns$ind
################################################################################
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

###############################################################################
# analyze rsa-states
tables.model <- states %>%
  dplyr::select(bn_id, r, probability, table.probs, table.support) %>%
  unnest(c(table.probs, table.support)) %>%
  group_by(bn_id) %>%
  pivot_wider(names_from = "table.support", values_from = "table.probs") %>%
  get_assertable_utterances()
  
theta = 0.80 
# when theta<=0.5 there are no empiric tables in independent contexts where only
# conditionals (and might) are assertable
utts.assertable = tables.model %>% 
  mutate(u_assertable = case_when(str_detect(utt, "might") ~ T, 
                                  prob >= theta ~ T, T ~ F)) %>% 
  filter(u_assertable)

utts.assertable %>% group_by(utt) %>% dplyr::count() %>% arrange(desc(n)) %>% 
  ungroup() %>% mutate(proportion = round(n / sum(n), 2))

# given only states where r = A || C (independence), given these states
# the conditional `if blue, green` is equally informative as `green`, 
# because whenever the conditional is assertable the literal is as well 
# (and vice versa) since P(g|b)=P(g) when r=A||C

# no literal assertable (-> no conjunction either, only 
# possibly a conditional and might)
bn_ids.only_might_or_if = utts.assertable %>% 
  mutate(might = utt %in% standardized.might) %>% 
  mutate(might_or_if = utt %in% standardized.might | utt %in% standardized.ifs) %>%
  group_by(bn_id) %>%
  mutate(s.might_or_if = sum(might_or_if), 
         s.might = sum(might), 
         n = n(), 
         only.might = s.might == n, 
         only.might_or_if = s.might_or_if == n) %>%
  filter(only.might_or_if) %>%
  pull(bn_id) %>% unique()

# for pe-task data: check utterance assertable in original PE-task tables
# tables.empiric <- data.observed %>%
#   dplyr::select(prolific_id, id, AC, `A-C`, `-AC`, `-A-C`) %>%
#   get_assertable_utterances()
# utts.assertable = tables.empiric %>% 
#   mutate(u_assertable = case_when(str_detect(utt, "might") ~ T, 
#                                   prob >= theta ~ T, T ~ F)) %>% 
#   filter(u_assertable)
# tbls.empiric.only_might_or_if <- utts.assertable %>% 
#   mutate(might = utt %in% standardized.might) %>% 
#   mutate(might_or_if = utt %in% standardized.might | utt %in% standardized.ifs) %>%
#   group_by(prolific_id, id) %>% 
#   mutate(s.might_or_if = sum(might_or_if), 
#          s.might = sum(might), 
#          n = n(), 
#          only.might = s.might == n, 
#          only.might_or_if = s.might_or_if == n) %>%
#   group_by(id, only.might_or_if) %>% dplyr::count() %>% arrange(desc(n)) %>% 
#   group_by(id) %>% mutate(proportion = n / sum(n)) %>% 
#   filter(str_detect(id, "independent")) %>% arrange(id, proportion) %>% 
#   filter(only.might_or_if)
# tbls.empiric.only_might_or_if
###############################################################################
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

save_data(weights %>% dplyr::select(id, bn_id, probs, prior_r, cn) %>% 
            add_column(n_rsa_states = params$n_forward_samples), 
          paste(result_dir, "weights_ci.rds", sep=FS))

# analyze weights
# check how likely states where only conditional or might are assertable are 
# for independent contexts
weights %>% 
  group_by(id, prior_r, cn) %>% 
  filter(bn_id %in% bn_ids.only_might_or_if & str_detect(id, "independent")) %>% 
  summarize(p = sum(probs), .groups = "drop_last") %>% arrange(prior_r, desc(p))

# relations for uninformative prior over relations
weights.relations = weights %>% filter(prior_r == "uninformative") %>% 
  group_by(id, r, probability) %>% 
  summarize(p = sum(probs), .groups = "drop_last") %>% 
  arrange(id, desc(p)) %>% mutate(p = round(p, 2)) %>% 
  filter(p != 0)
weights.relations %>% filter(id == "independent_ul")
weights.relations %>% filter(startsWith(id, "if1"))
weights.relations %>% filter(startsWith(id, "if2"))


weights.relations %>% 
  ggplot(aes(fill=probability)) + 
  geom_bar(aes(y=r, x=p), stat="identity") +
  facet_wrap(~id, scales="free")

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

  p <- data.joint %>% 
    filter(id == !!id) %>%
    filter(!(data=="empirical slider ratings"& prior_r=="informative")) %>% 
    mutate(data = case_when(data == "model states" ~ paste(prior_r, data), 
                            T ~ data)) %>% 
    ggplot(aes(x = cell)) +
    geom_bar(aes(y = ev_cell), stat = "identity") +
    facet_wrap(~data, scales = "free_x", 
               labeller=labeller(cell=cell_names), ncol=3) +
    labs(x = "event", y = "Expected value states/slider ratings", title = id)
  return(p)
})
plots.expected_states


N_obs = params$observations %>% group_by(id) %>% dplyr::count() 
N_rep = 100

map(slider_ratings$id %>% unique(), function(id){
  N = N_obs %>% filter(id == !!id) %>% pull(n)
  W <- weights %>% filter(id==!!id & prior_r == "informative") %>% rowid_to_column()
  
  # sample N_rep times same number of states as participants in this trial/context
  # and compute the expected value for each cell in each repetition
  evs = map_dfr(seq(1, N_rep), function(i){
    idx_states = sample(W$rowid, N, replace=T, prob=W$probs)
    states = W[idx_states,] %>% dplyr::select(rowid, id, AC, `A-C`, `-AC`, `-A-C`)
    expected.cells = states %>% 
      pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), names_to = "cell") %>% 
      group_by(cell) %>%
      summarize(ev_cell = mean(value), .groups = "drop_last") %>% 
      add_column(rep = i)
    return(expected.cells)
  })
  p <- slider_ratings %>% filter(id==!!id) %>% 
    ggplot() + 
    geom_boxplot(aes(x=cell, y=value)) +
    geom_boxplot(data=evs, aes(x=cell, y=ev_cell), color='firebrick') +
    stat_summary(aes(x=cell, y=value), fun=mean, geom="point", shape=4, size=5) +
    # geom_point(data = expected.states %>% filter(id==!!id), 
    #            aes(x=cell, y=ev_cell, color=prior_r), shape=4, stroke=2) +
    ggtitle(get_name_context(id))
  return(p)
})



