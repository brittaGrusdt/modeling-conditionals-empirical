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
library(latex2exp)
library(boot)
library(stringr)

source(here("R", "helpers-dirichlet-regression.R"))
source(here("R", "helpers-plotting.R"))

# Setup -------------------------------------------------------------------
theme_set(theme_clean(base_size = 20) + theme(legend.position = "top"))
prob_names <- c("blue"="P(b)", "green" = "P(g)")
# Data --------------------------------------------------------------------
active_config = "pathes"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "independent-contexts",
                   "dirichlet-regression-model", sep=FS)
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
brms_formula.ind <- y ~ pblue * pgreen + (1 + pblue * pgreen | subj) 
priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "mub"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "mug"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "munone"))

# Model independent contexts ----------------------------------------------
model_ind.pe_task = brm(
  data = df_ind.brms,
  family = "dirichlet",
  formula = brms_formula.ind,
  seed = 0710, 
  iter = 8000,
  chains = 4,
  cores = 4,
  warmup = 1000,
  prior = priors,
  control = list(adapt_delta = 0.9, max_treedepth=11),
  file = paste(target_dir, "model_pe_task_ind_dirichlet_maxreff.rds", sep = FS)
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
  
names_priors = c("H"="high", "L"="low", "U"="unc")
fn_pp_plots = function(trial_id){
  df.trial <- df.brms %>% filter(id == trial_id)
  N = nrow(df.trial)
  pb = names_priors[[str_sub(get_name_context(trial_id), -2, -2)]]
  pg = names_priors[[str_sub(get_name_context(trial_id), -1, -1)]]
  indices_b <- df_ind.brms$pblue == pb
  indices_g <- df_ind.brms$pgreen == pg
    
  y = df_ind.brms$y[indices_b & indices_g, ]
  y_vec <- matrix(y, ncol = N*4, byrow=T) %>% as.numeric()
  grp <- factor(c(rep("bg", N), rep("b¬g", N), rep("¬bg", N), rep("¬b¬g", N)), 
                levels = c("bg", "b¬g", "¬bg", "¬b¬g"), 
                labels = c("bg", "b¬g", "¬bg", "¬b¬g"))
  
  yrep = pp_samples[1:100, indices_b & indices_g, ]
  yrep_2d <- matrix(yrep, nrow=100, ncol=N*4)
  
  p <- ppc_dens_overlay_grouped(y=y_vec, yrep = yrep_2d, group = grp) +
    labs(title = get_name_context(trial_id)) +
    theme(panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(size=10),
          legend.position = "none") + 
    #by default ppc_dens_overlay uses fixed scales!
    facet_wrap("group", scales = "free") 
  
  fn <- paste("pp-tables-", trial_id, ".png", sep="")
  target_path <- paste(target_dir, FS, fn, sep="")
  ggsave(target_path, plot = p)
  return(p)
}
pp_plots <- map(df.brms %>% filter(relation=="independent") %>% pull(id) %>%
                  unique(), fn_pp_plots)


# bootstrapped data
get_mean_estimates = function(data, i){
  d2 <- data[i,] %>% as.matrix()
  return(colMeans(d2))
}
N = 1000
data_bootstrapped = group_map(
  df_ind.brms %>% group_by(pblue, pgreen), function(df.grp, grp){
    
    samples <- boot(df.grp$y %>% as_tibble(),
                    statistic = get_mean_estimates, R=N)
    bg = boot.ci(samples, index=c(1), type = "basic", conf=c(0.95))$basic[4:5]
    b = boot.ci(samples, index=c(2), type = "basic", conf=c(0.95))$basic[4:5]
    g = boot.ci(samples, index=c(3), type = "basic", conf=c(0.95))$basic[4:5]
    none = boot.ci(samples, index=c(4), type = "basic", conf=c(0.95))$basic[4:5]
    
    cis = rbind(bg, b, g, none)
    colnames(cis) <- c("lower", "upper")
    cis %>% as_tibble(rownames=NA) %>% rownames_to_column("world") %>% 
      add_column(pblue = grp$pblue, pgreen=grp$pgreen)
  }
) %>% bind_rows() %>%  add_column(data = "empiric") %>% 
  mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                        labels = c("bg", "b¬g", "¬bg", "¬b¬g")))

# plots with mean values
cols <- c("empiric" = "firebrick", "posterior predictive" = "gray")
pp_plots_means <- map(c("high", "low"), function(pgreen){
  map(c("high", "unc", "low"), function(pblue){
    if(!(pgreen == "high" && pblue == "low")){
      indices_blue <- df_ind.brms$pblue == pblue
      indices_green <- df_ind.brms$pgreen == pgreen
      indices <- indices_blue & indices_green
      
      pb <- switch(pblue, "high"="H", "unc"="U", "low"="L")
      pg <- switch(pgreen, "high"="H", "low"="L")
      tit <- paste("ind:", pb, pg, sep="")
        
      means_empiric = colMeans(df_ind.brms$y[indices, ]) %>% 
        enframe(name = "world", value = "mean_estimate")
      means_pp = colMeans(colMeans(pp_samples[, indices, ])) %>% 
        enframe(name = "world", value = "mean_estimate")
      
      df <- bind_rows(means_empiric %>% add_column(data = "empiric"),
                      means_pp %>% add_column(data = "posterior predictive")) %>% 
        mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                              labels = c("bg", "b¬g", "¬bg", "¬b¬g")))
      
      # compute highest posterior predictive density interval
      pp_means <- apply(pp_samples[,indices,], c(1,3), mean)
      pp_hdis <- apply(pp_means, c(2), tidybayes::hdi, .width=0.95 ) %>% 
       as_tibble() %>% add_column(limit=c("lower", "upper")) %>% 
        pivot_longer(cols = -limit, names_to = "world", values_to="val") %>% 
        pivot_wider(names_from = "limit", values_from = "val") %>% 
        add_column(data = "posterior predictive") %>% 
        mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                              labels = c("bg", "b¬g", "¬bg", "¬b¬g")))
      
      p <- df %>% 
        ggplot(aes(x = world, group = data, color=data, fill=data)) +
        geom_point(aes(y = mean_estimate)) + geom_line(aes(y = mean_estimate)) + 
        geom_ribbon(data = data_bootstrapped %>% 
                      filter(pblue == !!pblue & pgreen == !!pgreen), 
                    aes(ymin=lower, ymax=upper), alpha=0.5) +
        scale_color_manual(values = cols) +
        geom_ribbon(data = pp_hdis,aes(ymin=lower, ymax=upper), alpha=0.5) +
        labs(title = tit)
        
      ggsave(paste(target_dir, FS, "pp_", str_replace(tit, ":", "_"), ".png", sep=""), plot = p)
    return(p)
    }
  })
})

