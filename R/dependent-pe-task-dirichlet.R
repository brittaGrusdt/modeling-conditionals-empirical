library(here)
library(tidyverse)
library(tibble)
library(broom)
library(broom.mixed)
library(knitr)
library(ExpDataWrangling)
library(ModelUtils)
library(rwebppl)
library(bayesplot)
library(xtable)
library(ggpubr)
library(tidyselect)
library(bayesplot)
library(ggthemes)
library(tidybayes)
library(emmeans)
library(brms)

source(here("R", "helpers-dirichlet-regression.R"))
# Setup -------------------------------------------------------------------
p_cols = c("blue" = "blue4",
                  "green" = "forestgreen", 
                  "if_bg" = "hotpink1", 
                  "if_gb" = "sienna1" , 
                  "if_nbg" = "deeppink3", 
                  "if_ngb" = "orangered3", 
                  "AC" = "deepskyblue", 
                  "A-C" = "mediumblue", 
                  "-AC" = "darkorange", 
                  "-A-C" = "lightcoral")

dep_r_names <- c("if1" = "IF[1]", "if2" = "IF[2]")
theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
# Data --------------------------------------------------------------------
active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dep_dir = paste(here(params$dir_results), "dependent-contexts", sep=FS)
target_ind_dir = paste(here(params$dir_results), "independent-contexts", sep=FS)
if(!dir.exists(target_dep_dir)) dir.create(target_dep_dir, recursive = T)
if(!dir.exists(target_ind_dir)) dir.create(target_ind_dir, recursive = T)

data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task.smooth, slider) %>% 
  translate_standardized2model() 

data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task.smooth) %>% 
  pivot_wider(names_from = "utt.standardized", values_from = "pe_task.smooth") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`) %>% 
  get_controlled_factors() %>% dplyr::select(-relation_type) %>% 
  rename(blue = `blue falls`, green = `green falls`)


###########################################################################
######################  Dependent contexts  ###############################
###########################################################################
#  fit 4 slider responses dependent contexts  -----------------------------
df.brms <-  data.pe %>%
  dplyr::select(prolific_id, id, 
                AC, `A-C`, `-AC`, `-A-C`, 
                prior_blue, prior_green, relation) %>% 
  rename(bg = AC, b=`A-C`, g=`-AC`, none=`-A-C`, 
         pblue = prior_blue, pgreen = prior_green,
         subj = prolific_id) %>% 
  ungroup() %>% 
  mutate(s = (bg + b + g + none), 
         bg = bg/s, b = b/s, g = g/s, none = none/s)
df.brms$y <- with(df.brms, cbind(bg, b, g, none))
  
df_dep.brms <- df.brms %>% dplyr::select(subj, pblue, relation , y) %>% 
  mutate(relation = as.character(relation)) %>% 
  filter(relation != "independent") %>% 
  mutate(pblue = factor(pblue, levels = c("high", "unc", "uncl", "low")))
  
brms_formula.dep <- y ~ pblue * relation + (1 + pblue + relation | subj) 
# look at default priors
# get_prior(brms_formula.dep, data = df_dep.brms, family = "dirichlet")

# constrain priors for regression coefficients
priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b", dpar = "mub"),
            set_prior("normal(0, 1)", class = "b", dpar = "mug"),
            set_prior("normal(0, 1)", class = "b", dpar = "munone"))
            

model_dep.pe_task = brm(
  data = df_dep.brms,
  family = "dirichlet",
  formula = brms_formula.dep,
  seed = 0710, 
  iter = 4000,
  chains = 4,
  cores = 4,
  warmup = 1000,
  prior = priors,
  control = list(adapt_delta = 0.9, max_treedepth=11),
  file = paste(target_dep_dir, "model_pe_task_dep_dirichlet.rds", sep = FS)
)

posterior_draws <- tidy_draws(model_dep.pe_task)

# Analyze Dirichlet regression --------------------------------------------

# linear predictors mean estimate and quantiles
predictions <- map_dfr(c('Estimate', 'Q2.5', 'Q97.5'), function(val){
  model_params <- fixef(model_dep.pe_task)[, val]
  predicted_probs <- bind_rows(
    get_linear_predictors_dep('b', model_params), 
    get_linear_predictors_dep('g', model_params), 
    get_linear_predictors_dep('none', model_params)
  ) %>% 
    add_column(rowid = 1) %>% 
    pivot_longer(cols=c(-response_cat, -rowid), names_to = "predictor", 
                 values_to = "eta") %>% 
    transform_to_probs('bg') %>%
    pivot_wider(names_from = "response_cat", values_from = "p_hat") %>% 
    add_column(val = !!val)
  return(predicted_probs)
})

cat <- "bg"
predictions %>% group_by(predictor, val) %>% mutate(blue = bg + b) %>% 
  dplyr::select(predictor, val, !!cat) %>% 
  pivot_wider(names_from = "val", values_from = !!cat) %>% 
  add_column(response_cat = !!cat)

conds <- make_conditions(df_dep.brms, c("relation"))
p.cond_effects = conditional_effects(
  model_dep.pe_task, effects=c("pblue"), conditions = conds, 
  categorical = T, re_formula = NA, method = "fitted", robust = F
)
p.cond_effects

posterior_samples.probs = bind_rows(
  get_linear_predictors_samples_dep(posterior_draws, 'b'),
  get_linear_predictors_samples_dep(posterior_draws, 'g'),
  get_linear_predictors_samples_dep(posterior_draws, 'none')
) %>% pivot_longer(cols=c(-response_cat, -rowid), 
                   names_to = "predictor", values_to="eta") %>% 
  transform_to_probs('bg') %>% 
  pivot_wider(names_from = "response_cat", values_from = "p_hat") %>% 
  group_by(predictor) %>% 
  mutate(blue = bg + b) %>% 
  separate(predictor, into = c("relation", "prior_blue"), sep = "_")


# Plot posterior ----------------------------------------------------------
posterior_blue = posterior_samples.probs %>% 
  ggplot(aes(x=blue, fill=prior_blue)) + 
  stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
  facet_wrap(~relation, labeller = labeller(relation = dep_r_names, 
                                            .default = label_parsed)) +
  labs(x = "P(blue)", y = "")
posterior_blue
ggsave(paste(target_dep_dir, "posterior_blue.png", sep=FS), plot = posterior_blue)

posterior_green = posterior_samples.probs %>% 
  mutate(green = bg + g) %>% 
  ggplot(aes(x=green, fill=prior_blue)) + 
  stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
  facet_wrap(~relation, labeller = labeller(relation = dep_r_names, 
                                            .default = label_parsed)) +
  labs(x = "P(green)", y = "")
posterior_green
ggsave(paste(target_dep_dir, "posterior_green.png", sep=FS), plot = posterior_green)

posterior_if_bg = posterior_samples.probs %>% 
  mutate(if_bg = bg/blue) %>% 
  ggplot(aes(x=if_bg, fill=prior_blue)) + 
  stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
  facet_wrap(~relation, labeller = labeller(relation = dep_r_names, 
                                            .default = label_parsed)) +
  labs(x = "P(g|b)", y = "")
posterior_if_bg
ggsave(paste(target_dep_dir, "posterior_if_bg.png", sep=FS), plot = posterior_if_bg)

posterior_if_nbg = posterior_samples.probs %>% 
  mutate(if_nbg = g/(g+none)) %>% 
  ggplot(aes(x=if_nbg, fill=prior_blue)) + 
  stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
  facet_wrap(~relation, labeller = labeller(relation = dep_r_names, 
                                            .default = label_parsed)) +
  labs(x = "P(g|¬b)", y = "")
posterior_if_nbg
ggsave(paste(target_dep_dir, "posterior_if_nbg.png", sep=FS), plot = posterior_if_nbg)


# posterior probabilities hypotheses
# P(blue) estimated given prior blue: high>unc>uncl>low
cat <- "bg"
posterior_samples.probs %>% 
  dplyr::select(rowid, relation, prior_blue, !!cat) %>% 
  group_by(rowid, relation) %>% 
  pivot_wider(names_from = "prior_blue", values_from = !!cat) %>% 
  group_by(relation) %>% 
  summarize(high_unc = mean(high > unc), 
            unc_uncl = mean(unc > uncl),
            unc_low = mean(unc > low),
            uncl_low = mean(uncl > low)) %>% 
  add_column(response = !!cat)



# Posterior predictive checks ---------------------------------------------
pp_samples <- posterior_predict(model_dep.pe_task)

cat <- "b"
ppc_dens_overlay_grouped(y=df_dep.brms$y[,cat], 
                         yrep = pp_samples[1:100,,cat], 
                         group = interaction(df_dep.brms$relation, 
                                             df_dep.brms$pblue))








######## DIRICHLET REGRESSION TESTS #######################################

# Test Dirichlet regression -----------------------------------------------
y1 <- rdirichlet(1000, c(8, 1, 1, 8))
colnames(y1) <- c("bg", "b", "g", "none")
y2 = rdirichlet(1000, c(1, 8, 8, 1))
colnames(y2) <- c("bg", "b", "g", "none")

df.sim <- bind_rows(as_tibble(y1) %>% add_column(x = "A"),
                    as_tibble(y2) %>% add_column(x = "B"))
df.sim$y <- with(df.sim, cbind(bg, b, g, none))

df.sim.mean = df.sim %>% group_by(x) %>% 
  summarize(mean.bg = mean(bg), 
            mean.b = mean(b),
            mean.g = mean(g),
            mean.none = mean(none))
df.sim.mean

model.dir_test <- brm(
  data = df.sim,
  family = "dirichlet",
  formula = y ~ 1 + x,
  #seed = 0710, 
  iter = 2000,
  chains = 4,
  cores = 4,
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth=11),
  file = "model_test_dirichlet.rds"
)
model_params <- fixef(model.dir_test)[, 'Estimate']

get_test_predictor = function(predicted_cat, model_params){
  mu = paste('mu', predicted_cat, sep='')
  pars = model_params[startsWith(names(model_params), mu)] 
  df = tibble(
    prediction = predicted_cat,
    A = pars[[paste(mu, '_Intercept', sep='')]],
    B = A + pars[[paste(mu, '_xB', sep = '')]]
  )
  return(df)
}

bind_rows(
  get_test_predictor('b', model_params), 
  get_test_predictor('g', model_params), 
  get_test_predictor('none', model_params)
) %>% 
  pivot_longer(cols=c(-prediction), names_to = "param") %>% 
  mutate(p_eta = exp(value)) %>% 
  group_by(param) %>% 
  mutate(denominator = 1 + sum(p_eta),
         p_refcat = 1 / denominator, 
         p_cat = p_refcat * p_eta) %>% filter(param == "A")
