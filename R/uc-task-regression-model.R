library(brms)
library(here)
library(tidyverse)
library(ExpDataWrangling)
library(bayesplot)
library(boot)
library(ggthemes)
library(latex2exp)

theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))

dep_r_names <- c("if1" = "IF[1]", "if2" = "IF[2]")

# Data --------------------------------------------------------------------
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

plot_dir = here(params$dir_data, "figs")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = T)

path_cleaned_data = here("data", "cleaned-data.csv")
data = get_controlled_factors(read_csv(path_cleaned_data)) %>% 
  mutate(utterance = utt.standardized) %>% 
  chunk_utterances() %>%
  rename(utt_type = utterance, utterance = utt.standardized) %>% 
  mutate(utt_type = as.character(utt_type))

data.uc = data %>% 
  filter(human_exp2 == 1) %>% 
  dplyr::select(prolific_id, id, relation, relation_type, prior_blue, prior_green,
                utt_type, uc_task,pe_task.smooth, pe_task) %>% 
  mutate(conditional = case_when(utt_type != "conditional" ~ 0, T ~ 1),
         # prior_blue_agg = case_when(as.character(prior_blue) %in% c("unc", "uncl") ~ "unc",
         #                            T ~ "confident"),
         prior_blue_agg = case_when(as.character(prior_blue) == "unc" ~ "unc", 
                                    as.character(prior_blue) == "uncl" ~ "uncl",
                                    T ~ "confident"),
         prior_blue_agg = factor(prior_blue_agg, levels = c("unc", "uncl", 
                                                            "confident"),
                                 labels = c("unc", "uncl", "confident")),
         prior_blue = factor(prior_blue,
                             levels = c("low", "uncl", "unc", "high"),
                             labels = c("L", "U^-", "U", "H"))
         ) %>% 
  rename(subject_id = prolific_id)

df <- data.uc %>% 
  mutate(conditional = as.logical(conditional), 
         prior_blue_agg = case_when(as.character(prior_blue) %in% c("U", "U^-") ~ "uncertain",
                                    as.character(prior_blue) == "H" ~ "high", 
                                    as.character(prior_blue) == "L" ~ "low", 
                                    T ~ as.character(prior_blue)))
model_bernoulli = brm(
  data = df,
  family = "bernoulli",
  formula = conditional ~ 1 + prior_blue_agg * relation + 
    (1 + prior_blue_agg * relation |subject_id),
  seed = 0710, 
  warmup = 1000,
  iter = 8000,
  control = list(adapt_delta = 0.95, max_treedepth=12),
  file = here("results", "model-uc-task-bernoulli.rds"),
  chains = 4,
  cores = 4
)
pp_samples <- posterior_predict(model_bernoulli)
ppc_dens_overlay_grouped(y=df$conditional %>% as.numeric(),
                         yrep = pp_samples[1:500,],
                         group = interaction(df$relation, 
                                             df$prior_blue_agg))

# if1 uncertain vs. confident
if1_high <- "Intercept + relationif1"
if1_unc <- "Intercept + relationif1 + prior_blue_agguncertain + prior_blue_agguncertain:relationif1"
if1_low <- "Intercept + relationif1 + prior_blue_agglow + prior_blue_agglow:relationif1"
h1 <- paste(if1_unc, ">", if1_low)
hypothesis(model_bernoulli, h1)
h2 <- paste(if1_unc, ">", if1_high)
hypothesis(model_bernoulli, h2)

# if2 uncertain vs. confident
if2_high <- "Intercept"
if2_unc <- "Intercept + prior_blue_agguncertain"
if2_low <- "Intercept + prior_blue_agglow"
hypothesis(model_bernoulli, paste(if2_unc, ">", if2_high))
hypothesis(model_bernoulli, paste(if2_unc, ">", if2_low))

# independent uncertain vs. confident
ind_high <- "inv_logit_scaled(Intercept + relationindependent)"
ind_low <- "inv_logit_scaled(Intercept + relationindependent + prior_blue_agglow + prior_blue_agglow:relationindependent)"
ind_unc <- "inv_logit_scaled(Intercept + relationindependent + prior_blue_agguncertain + prior_blue_agguncertain:relationindependent)"

# independent == 0
hypothesis(model_bernoulli, paste(ind_unc, "= 0"))
hypothesis(model_bernoulli, paste(ind_high, "= 0"))
hypothesis(model_bernoulli, paste(ind_low, "= 0"))



