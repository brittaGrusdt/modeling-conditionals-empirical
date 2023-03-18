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
config_cns = "fine_grained_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "semi_informative"
params <- prepare_data_for_wppl(config_cns, config_weights_relations, 
                                extra_packages = extra_packages)

slider_ratings = params$observations %>% ungroup() %>% 
  dplyr::select(id, prolific_id, AC, `A-C`, `-AC`, `-A-C`) %>% 
  pivot_longer(cols = c("AC", "A-C", "-AC", "-A-C"), 
               names_to = "cell", values_to = "value")

###############################################################################
# analyze RSA-states
theta = 0.894
utts.assertable <- params$prior_samples %>%
  dplyr::select(bn_id, r, probability, table.probs, table.support) %>%
  unnest(c(table.probs, table.support)) %>%
  group_by(bn_id) %>%
  pivot_wider(names_from = "table.support", values_from = "table.probs") %>%
  get_assertable_utterances(theta)
  
# how often is each utterance assertable (in how many states)?
utts.assertable %>% group_by(utt) %>% dplyr::count() %>% arrange(desc(n)) %>% 
  ungroup() %>% mutate(proportion = round(n / params$nb_rsa_states , 2))

# given only states where r = A || C (independence), given these states
# the conditional `if blue, green` is equally informative as `green`, 
# because whenever the conditional is assertable the literal is as well 
# (and vice versa) since P(g|b)=P(g) when r=A||C

# no literal assertable (-> no conjunction either, only 
# possibly a conditional and might)
bn_ids.only_might_and_if = utts.assertable %>% 
  mutate(might = utt %in% standardized.might) %>% 
  mutate(might_or_if = utt %in% standardized.might | utt %in% standardized.ifs) %>%
  group_by(bn_id) %>%
  mutate(s.might_or_if = sum(might_or_if), 
         s.might = sum(might), 
         n = n(), 
         only.might = s.might == n, 
         only.might_or_if = s.might_or_if == n) %>%
  filter(only.might_or_if & !only.might) %>%
  pull(bn_id) %>% unique()

# analyze assertable utterances PE-task data
# for PE-task data: check utterance assertable in original PE-task tables
utts.assertable.emp <- params$observations %>%
  dplyr::select(prolific_id, id, AC, `A-C`, `-AC`, `-A-C`, utterance) %>%
  rename(uc_task = utterance) %>%
  get_assertable_utterances(theta)
tbls.empiric.only_might_or_if <- utts.assertable.emp %>%
  mutate(might = utt %in% standardized.might) %>%
  mutate(might_or_if = utt %in% standardized.might | utt %in% standardized.ifs) %>%
  group_by(prolific_id, id) %>%
  mutate(s.might_or_if = sum(might_or_if),
         s.might = sum(might),
         n = n(),
         only.might = s.might == n,
         only.might_or_if = s.might_or_if == n) %>%
  group_by(id, only.might_or_if, only.might) %>% dplyr::count() %>% arrange(desc(n)) %>%
  group_by(id) %>% mutate(proportion = n / sum(n)) %>%
  filter(str_detect(id, "independent")) %>% arrange(id, proportion) %>%
  filter(only.might_or_if) %>% dplyr::select(-only.might_or_if)
tbls.empiric.only_might_or_if
###############################################################################


# compute weights ---------------------------------------------------------
weights = map(c("informative", "semi_informative"), function(prior_r){
  
  params <- prepare_data_for_wppl(config_cns, prior_r, 
                                  extra_packages = extra_packages)
  params$verbose <- FALSE
  weights_ci = run_webppl("webppl-model/weights-contexts.wppl", params) %>% 
    bind_rows(.id = "id") %>% group_by(id) %>% arrange(desc(probs)) %>% 
    unnest(c(support)) %>% 
    unnest(c(table.probs, table.support)) %>% 
    pivot_wider(names_from = table.support, values_from = table.probs) %>% 
    group_by(id) %>% mutate(cdf = cumsum(probs)) %>% 
    add_column(prior_r = prior_r)
  
  save_data(weights_ci %>% add_column(n_rsa_states = params$nb_rsa_states),
            paste(params$config_dir, "weights_ci.rds", sep = FS))
  return(weights_ci)
}) %>% bind_rows() 


# analyze weights ---------------------------------------------------------
# check how likely states where only conditional or might are assertable are 
# for independent contexts
# (should not be the case at all if prior over relations =1 for A||C)
weights %>% 
  group_by(id, prior_r, cn) %>% 
  filter(bn_id %in% bn_ids.only_might_and_if & str_detect(id, "independent")) %>% 
  summarize(p = sum(probs), .groups = "drop_last") %>% arrange(prior_r, desc(p))

# relations for uninformative prior over relations
weights.relations = weights %>% filter(prior_r == "semi_informative") %>% 
  group_by(id, r, probability) %>% 
  summarize(p = sum(probs), .groups = "drop_last") %>% 
  arrange(id, desc(p)) %>% mutate(p = round(p, 2)) %>% 
  filter(p != 0)
weights.relations %>% filter(startsWith(id, "if1"))
weights.relations %>% filter(startsWith(id, "if2"))


weights.relations %>% filter(str_detect(id, "if")) %>% 
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
  expected.slider_ratings %>% add_column(prior_r = "semi_informative")
) %>% 
  mutate(cell = factor(cell, levels = c("AC", "A-C", "-AC", "-A-C")))



# some plots --------------------------------------------------------------
# plot expected ratings / model states
plots.expected_states = map(data.joint$id %>% unique, function(id){
  p1 <- data.joint %>%
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

  p2 <- data.joint %>% 
    filter(id == !!id) %>%
    filter(!(data=="empirical slider ratings"& prior_r=="informative")) %>% 
    mutate(data = case_when(data == "model states" ~ paste(prior_r, data), 
                            T ~ data)) %>% 
    ggplot(aes(x = cell)) +
    geom_bar(aes(y = ev_cell), stat = "identity") +
    facet_wrap(~data, scales = "free_x", 
               labeller=labeller(cell=cell_names), ncol=3) +
    labs(x = "event", y = "Expected value states/slider ratings", title = id)
  return(p2)
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



