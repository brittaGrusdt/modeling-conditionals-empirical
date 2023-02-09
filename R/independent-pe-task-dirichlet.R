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

theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
# Data --------------------------------------------------------------------
active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "independent-contexts", sep=FS)
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)

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

###########################################################################
#####################  Independent contexts  ##############################
###########################################################################
#  fit 4 slider responses independent contexts  ----------------------------

# Data --------------------------------------------------------------------
df_ind.brms <-  df.brms %>% dplyr::select(subj, pblue, pgreen, relation, y) %>% 
  mutate(relation = as.character(relation), 
         pblue = as.character(pblue),
         pgreen = as.character(pgreen)) %>% 
  filter(relation == "independent")
brms_formula.ind <- y ~ pblue * pgreen + (1 + pblue + pgreen | subj) 
priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b", dpar = "mub"),
            set_prior("normal(0, 1)", class = "b", dpar = "mug"),
            set_prior("normal(0, 1)", class = "b", dpar = "munone"))

# Model independent contexts ----------------------------------------------
model_ind.pe_task = brm(
  data = df_ind.brms,
  family = "dirichlet",
  formula = brms_formula.ind,
  seed = 0710, 
  iter = 4000,
  chains = 4,
  cores = 4,
  warmup = 1000,
  prior = priors,
  control = list(adapt_delta = 0.9, max_treedepth=11),
  file = paste(target_dir, "model_pe_task_ind_dirichlet.rds", sep = FS)
)

posterior_draws <- tidy_draws(model_ind.pe_task)

# Analyze Dirichlet regression --------------------------------------------
# linear predictors mean estimate
model_params <- fixef(model_ind.pe_task)[, 'Estimate']
predictions <- bind_rows(
  get_linear_predictors_ind('b', model_params), 
  get_linear_predictors_ind('g', model_params), 
  get_linear_predictors_ind('none', model_params)
) %>% 
  add_column(rowid = 1) %>% 
  pivot_longer(cols=c(-response_cat, -rowid), names_to = "predictor", 
               values_to = "eta") %>% 
  transform_to_probs('bg') %>%
  pivot_wider(names_from = "response_cat", values_from = "p_hat")

conds <- make_conditions(df_ind.brms, c("pgreen"))
p.cond_effects = conditional_effects(
  model_ind.pe_task, effects=c("pblue"), conditions = conds, 
  categorical = T, re_formula = NA, method = "fitted", robust = F
)
p.cond_effects

posterior_samples.probs = bind_rows(
  get_linear_predictors_samples_ind(posterior_draws, 'b'),
  get_linear_predictors_samples_ind(posterior_draws, 'g'),
  get_linear_predictors_samples_ind(posterior_draws, 'none')
) %>% pivot_longer(cols=c(-response_cat, -rowid), 
                   names_to = "predictor", values_to="eta") %>% 
  transform_to_probs('bg') %>% 
  pivot_wider(names_from = "response_cat", values_from = "p_hat") %>% 
  group_by(predictor) %>% 
  mutate(blue = bg + b, green = bg + g) %>% 
  separate(predictor, into = c("relation", "prior_blue_green"), sep = "_")


# Plot posterior ----------------------------------------------------------
plot_posterior_prob = function(posterior_samples.probs, x, lab_x){
  plot.posterior = posterior_samples.probs %>% 
    ggplot(aes(x=get(x), fill=prior_blue_green)) + 
    stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
    labs(x = lab_x, y = "") 
  return(plot.posterior)
}

posterior_blue = posterior_samples.probs %>% plot_posterior_prob("blue", "P(blue)")
ggsave(paste(target_dir, "posterior_blue.png", sep=FS), plot = posterior_blue)

posterior_green = posterior_samples.probs %>% plot_posterior_prob("green", "P(green)")
ggsave(paste(target_dir, "posterior_green.png", sep=FS), plot = posterior_green)

posterior_bg = posterior_samples.probs %>% 
  plot_posterior_prob("bg", "P(b,g)")
ggsave(paste(target_dir, "posterior_bg.png", sep=FS), plot = posterior_bg)

posterior_if_bg = posterior_samples.probs %>% mutate(if_bg = bg/(bg+b)) %>% 
  plot_posterior_prob("if_bg", "P(g|b)")
ggsave(paste(target_dir, "posterior_if_bg.png", sep=FS), plot = posterior_if_bg)

posterior_if_nbg = posterior_samples.probs %>% mutate(if_nbg = g/(g+none)) %>% 
  plot_posterior_prob("if_nbg", "P(g|¬b)")
ggsave(paste(target_dir, "posterior_if_nbg.png", sep=FS), plot = posterior_if_nbg)


posterior_if_gb = posterior_samples.probs %>%  mutate(if_gb = bg/(bg+g)) %>% 
  plot_posterior_prob("if_gb", "P(b|g)")
ggsave(paste(target_dir, "posterior_if_gb.png", sep=FS), plot = posterior_if_gb)

posterior_if_ngb = posterior_samples.probs %>% mutate(if_ngb = b/(b+none)) %>% 
  plot_posterior_prob("if_ngb", "P(b|¬g)")
ggsave(paste(target_dir, "posterior_if_ngb.png", sep=FS), plot = posterior_if_ngb)

# plot 4 probability sliders at once
posterior_cells = plot_posterior_prob(
  posterior_samples.probs %>% 
    dplyr::select(-blue, -green, -relation) %>% 
    pivot_longer(cols=c(bg, b, g, none), names_to = "world", values_to = "p_hat") %>% 
    mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                          labels = c("bg", "b¬g", "¬bg", "¬b¬g")),
           prior_blue_green=factor(prior_blue_green, 
                                   levels = c("ll", "ul", "hl", "uh", "hh"), 
                                   labels = c("LL", "UL", "HL", "UH", "HH")
                                   )),
  "p_hat", "Estimated probability"
) + facet_wrap(~world, ncol = 2, labeller = labeller(relation = label_parsed)) +
  theme(panel.spacing = unit(2, "lines"))
ggsave(paste(target_dir, "posterior_cells.png", sep=FS), 
       plot = posterior_cells, width = 12)



# posterior probabilities hypotheses --------------------------------------
# P(blue) estimated given prior blue: high>unc>uncl>low

cat <- "blue"
posterior_samples.probs %>% 
  dplyr::select(rowid, prior_blue_green, !!cat) %>% 
  group_by(rowid, prior_blue_green) %>% 
  pivot_wider(names_from = "prior_blue_green", values_from = !!cat) %>% 
  ungroup() %>% 
  summarize(ll_ul = mean(ll < ul), 
            ul_hh = mean(ul < hh),
            hl_uh = mean(hl > uh), 
            ul_uh = mean(ul < uh), 
            hh_hl = mean(hh < hl)
            ) %>% 
  add_column(response = !!cat)


cat <- "green"
  posterior_samples.probs %>% 
    dplyr::select(rowid, prior_blue_green, !!cat) %>% 
    group_by(rowid, prior_blue_green) %>% 
    pivot_wider(names_from = "prior_blue_green", values_from = !!cat) %>% 
    ungroup() %>% 
    summarize(ll_ul = mean(ll < ul), 
              ll_hl = mean(ll < hl),
              ul_hl = mean(ul < hl)) %>% 
    add_column(response = !!cat)
  


# Posterior predictive checks ---------------------------------------------
pp_samples <- posterior_predict(model_ind.pe_task)

pp_plots <- map(c("bg", "b", "g", "none"), function(cat){
  map(c("high", "unc", "low"), function(pblue){
    indices <- df_ind.brms$pblue == pblue
    tit <- switch(cat, "bg" = "bg", "b" = "b¬g", "g" = "¬bg", "none" = "¬b¬g")
    xlab <- switch(pblue, "high"="H", "unc"="U", "low"="L")
    p <- ppc_dens_overlay_grouped(
      y=df_ind.brms$y[indices, cat], 
      yrep = pp_samples[1:100, indices, cat], 
      group = df_ind.brms$pgreen[indices]
      ) +
      labs(x = paste("prior_blue:", xlab), title = tit) +
      theme(legend.position = "none")
    ggsave(paste(target_dir, FS, "pp_", cat, "_pblue_", pblue, 
                 ".png", sep=""), plot = p, width = 5, height=3)
    return(p)
  })
})













