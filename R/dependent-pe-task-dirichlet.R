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

source(here("R", "helpers-dirichlet-regression.R"))
source(here("R", "helpers-data-models.R"))
source(here("R", "helpers-plotting.R"))
# Setup -------------------------------------------------------------------
theme_set(theme_clean(base_size = 26) + theme(legend.position = "top"))

# Data --------------------------------------------------------------------
active_config = "pathes"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "dependent-contexts",
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
  
brms_formula.dep <- y ~ pblue * relation + (1 + pblue * relation | subj) 
# look at default priors
# get_prior(brms_formula.dep, data = df_dep.brms, family = "dirichlet")

# constrain priors for regression coefficients
priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "mub"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "mug"),
            set_prior("normal(0, 2.5)", class = "b", dpar = "munone"))
            

model_dep.pe_task = brm(
  data = df_dep.brms,
  family = "dirichlet",
  formula = brms_formula.dep,
  seed = 0710, 
  iter = 8000,
  chains = 4,
  cores = 4,
  warmup = 1000,
  prior = priors,
  control = list(adapt_delta = 0.9, max_treedepth=11),
  file = paste(target_dir, "model_pe_task_dep_dirichlet_maxreff.rds", sep = FS)
)

posterior_draws <- tidy_draws(model_dep.pe_task)

# Analyze Dirichlet regression --------------------------------------------

# linear predictors estimates
model_params <- fixef(model_dep.pe_task)[, 'Estimate']
predictions <- bind_rows(
  get_linear_predictors_dep('b', model_params), 
  get_linear_predictors_dep('g', model_params), 
  get_linear_predictors_dep('none', model_params)
) %>% 
  add_column(rowid = 1) %>% 
  pivot_longer(cols=c(-response_cat, -rowid), names_to = "predictor", 
               values_to = "eta") %>% 
  transform_to_probs('bg') %>%
  pivot_wider(names_from = "response_cat", values_from = "p_hat")

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
) %>% pivot_longer(cols=c(-response_cat, -rowid, -.chain, -.iteration), 
                   names_to = "predictor", values_to="eta") %>% 
  transform_to_probs('bg') %>% 
  pivot_wider(names_from = "response_cat", values_from = "p_hat") %>% 
  group_by(predictor) %>% 
  mutate(blue = bg + b) %>% 
  separate(predictor, into = c("relation", "prior_blue"), sep = "_") %>% 
  mutate(prior_blue = factor(prior_blue, levels = c("low", "uncl", "unc", "high")))


# Plot posterior ----------------------------------------------------------
plot_posterior_prob = function(posterior_samples.probs, x, lab_x){
  plot.posterior = posterior_samples.probs %>% 
    mutate(relation = factor(relation, levels = c("if1", "if2"),
                             labels = c(parse(text = expression("if"[1])),
                                        parse(text = expression("if"[2])))),
           prior_blue = factor(prior_blue,
                               levels = c("low", "uncl", "unc", "high"),
                               labels = c("L", "U^-", "U", "H"))) %>% 
    ggplot(aes(x=get(x), fill=prior_blue)) + 
    stat_halfeye(.width = c(0.95), point_interval = "mean_hdi") +
    facet_wrap(~relation, labeller = label_parsed) +
    scale_fill_discrete(labels = unname(TeX(c("L", "$U^-$", "U", "H")))) + 
    labs(x = lab_x, y = "") 
  return(plot.posterior)
}

# plot all 4 probability estimates in one plot:
posterior_cells = plot_posterior_prob(
  posterior_samples.probs %>% pivot_longer(cols=c(bg, b, g, none), 
                                           names_to = "world",
                                           values_to = "p_hat") %>% 
    mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                          labels = c("bg", "b¬g", "¬bg", "¬b¬g"))),
  "p_hat", "Estimated probability"
) + facet_grid(relation~world, labeller = labeller(relation = label_parsed)) +
  theme(panel.spacing = unit(2, "lines"))
posterior_cells
ggsave(paste(target_dir, "posterior_cells.png", sep=FS), 
       plot = posterior_cells, width = 12)


posterior_blue = plot_posterior_prob(posterior_samples.probs, "blue", "P(b)")
ggsave(paste(target_dir, "posterior_blue.png", sep=FS), plot = posterior_blue)

posterior_green =  plot_posterior_prob(posterior_samples.probs %>%
                                         mutate(green = bg + g),
                                       "green", "P(g)")
ggsave(paste(target_dir, "posterior_green.png", sep=FS), plot = posterior_green)

posterior_if_bg =  plot_posterior_prob(posterior_samples.probs %>% 
                                         mutate(if_bg = bg/blue), 
                                       "if_bg", "P(g|b)")
ggsave(paste(target_dir, "posterior_if_bg.png", sep=FS), plot = posterior_if_bg)

posterior_if_nbg =  plot_posterior_prob(posterior_samples.probs %>% 
                                          mutate(if_nbg = g/(g+none)),
                                        "if_nbg", "P(g|¬b)")
  ggsave(paste(target_dir, "posterior_if_nbg.png", sep=FS), plot = posterior_if_nbg)


# Hypotheses --------------------------------------------------------------
# P(bg) estimated given prior blue: high>unc>uncl>low
posterior_samples.probs %>% 
  dplyr::select(rowid, relation, prior_blue, "bg") %>% 
  group_by(rowid, relation) %>% 
  pivot_wider(names_from = "prior_blue", values_from = "bg") %>% 
  group_by(relation) %>% 
  summarize(high_unc = mean(high > unc),
            high_uncl = mean(high > uncl),
            unc_uncl = mean(unc > uncl),
            unc_low = mean(unc > low),
            uncl_low = mean(uncl > low)) %>% 
  add_column(response = "bg")


posterior_samples.probs %>% 
  dplyr::select(rowid, relation, prior_blue, "none") %>% 
  group_by(rowid, relation) %>% 
  pivot_wider(names_from = "prior_blue", values_from = "none") %>% 
  group_by(relation) %>% 
  summarize(high_unc = mean(high < unc), 
            high_uncl = mean(high < uncl), 
            unc_uncl = mean(unc < uncl),
            unc_low = mean(unc < low),
            uncl_low = mean(uncl < low)) %>% 
  add_column(response = "none")


# P(g|¬b)
posterior_samples.probs %>% 
  mutate(if_nbg = g/(g+none)) %>% 
  dplyr::select(rowid, relation, prior_blue, if_nbg) %>% 
  pivot_wider(values_from = if_nbg, names_from =relation) %>% 
  mutate(if2_larger = if2 > if1) %>% 
  group_by(prior_blue) %>% 
  summarize(p_h = sum(if2_larger)/n(), .groups = "drop_last")
  


# Posterior predictive checks ---------------------------------------------
pp_samples <- posterior_predict(model_dep.pe_task)
pp_plots <- map(df.brms %>% filter(relation != "independent") %>%
                  pull(id) %>% unique(), 
                function(trial_id){
  df.trial <- df.brms %>% filter(id == trial_id)
  tit_pblue <- get_name_context(trial_id)
  indices1 <- df_dep.brms$pblue == eval(tit_pblue[2])
  indices2 <- df_dep.brms$relation == substr(trial_id, 1, 3)
  N = nrow(df.trial)
  y = df_dep.brms$y[indices1 & indices2, ]
  y_vec <- matrix(y, ncol = N*4, byrow=T) %>% as.numeric()
  grp <- factor(c(rep("bg", N), rep("b¬g", N), rep("¬bg", N), rep("¬b¬g", N)), 
                levels = c("bg", "b¬g", "¬bg", "¬b¬g"), 
                labels = c("bg", "b¬g", "¬bg", "¬b¬g"))
  
  yrep = pp_samples[1:100, indices1 & indices2, ]
  yrep_2d <- matrix(yrep, nrow=100, ncol=N*4)

  p <- ppc_dens_overlay_grouped(y=y_vec, yrep = yrep_2d, group = grp) +
    labs(title = parse(text=tit_pblue[1])) +
    theme(panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(size=10),
          legend.position = "none") + 
    #by default ppc_dens_overlay uses fixed scales!
    facet_wrap("group", scales = "free") 
  
  fn <- paste(target_dir, FS, paste("pp-tables-", trial_id, ".png", sep=""), sep="")
  ggsave(fn, plot = p, width = 5, height=5)
  return(p)
})

# bootstrapped data
get_mean_estimates = function(data, i){
  d2 <- data[i,] %>% as.matrix()
  return(colMeans(d2))
}
N = 1000
data_bootstrapped = group_map(
  df_dep.brms %>% group_by(pblue, relation), function(df.grp, grp){
    
    samples <- boot(df.grp$y %>% as_tibble(),
                    statistic = get_mean_estimates, R=N)
    bg = boot.ci(samples, index=c(1), type = "basic", conf=c(0.95))$basic[4:5]
    b = boot.ci(samples, index=c(2), type = "basic", conf=c(0.95))$basic[4:5]
    g = boot.ci(samples, index=c(3), type = "basic", conf=c(0.95))$basic[4:5]
    none = boot.ci(samples, index=c(4), type = "basic", conf=c(0.95))$basic[4:5]
    
    cis = rbind(bg, b, g, none)
    colnames(cis) <- c("lower", "upper")
    cis %>% as_tibble(rownames=NA) %>% rownames_to_column("world") %>% 
      add_column(pblue = grp$pblue, relation = grp$relation)
  }) %>% bind_rows() %>%  add_column(data = "empiric") %>% 
  mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                        labels = c("bg", "b¬g", "¬bg", "¬b¬g")))

# plots with mean values
cols <- c("empiric" = "firebrick", "posterior predictive" = "gray")
pp_plots_means <- map(c("if1", "if2"), function(relation){
  map(c("high", "unc", "uncl", "low"), function(pblue){
    indices_rel <- df_dep.brms$relation == relation
    indices_blue <- df_dep.brms$pblue == pblue
    indices <- indices_blue & indices_rel
    cond <- paste(relation, pblue, sep="_")
    tit <- switch(cond, "if1_high"=expression("if"[1]*":HI"), 
                        "if1_unc"=expression(paste(`if`[1], ":UI")),
                        "if1_uncl"=expression("if"[1]*":U"^-{}*"I"),
                        "if1_low"=expression("if"[1]*":LI"),
                        "if2_high"=expression("if"[2]*":HL"), 
                        "if2_unc"=expression("if"[2]*":UL"),
                        "if2_uncl"=expression("if"[2]*":U"^-{}*"L"),
                        "if2_low"=expression("if"[2]*":LL"))
    
    means_empiric = colMeans(df_dep.brms$y[indices, ]) %>% 
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
                    filter(pblue == !!pblue & relation == !!relation), 
                  aes(ymin=lower, ymax=upper), alpha=0.5) +
      scale_color_manual(values = cols) +
      geom_ribbon(data = pp_hdis,aes(ymin=lower, ymax=upper), alpha=0.5) +
      ggtitle(tit)
    
      ggsave(paste(target_dir, FS, "pp-evs-", cond, ".png", sep=""), plot = p)
    return(p)
  })
})

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
