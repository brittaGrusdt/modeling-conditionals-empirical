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
library(ggthemes)

source(here("R", "helpers-data-models.R"))

# for plots
theme_set(theme_clean(base_size = 26) + 
            theme(legend.position = "top", 
                  axis.text = element_text(size = 26),
                  axis.title = element_text(size = 26)))
prob_names <- c("blue"="P(b)", "green" = "P(g)")
trial_names <- c("independent_hh" = "ind:HH", 
                 "independent_ll" = "ind:LL", 
                 "independent_uh" = "ind:UH", 
                 "independent_ul" = "ind:UL", 
                 "independent_hl" = "ind:HL")

# Data --------------------------------------------------------------------
active_config = "pathes"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
params <- config::get()

target_dir = paste(here(params$dir_results), "independent-contexts", 
                   "zoib-model", sep=FS)
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, 
                pe_task.smooth, pe_task, slider) %>% 
  translate_standardized2model() 

# use unsmoothed slider ratings from PE-task! 
# (data model with zero-one inflated beta!)
data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task) %>% 
  pivot_wider(names_from = "utt.standardized", 
              values_from = "pe_task") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`) %>% 
  get_controlled_factors() %>% dplyr::select(-relation_type) %>% 
  rename(blue = `blue falls`, green = `green falls`) %>% 
  mutate(blue_zero_one = (blue == 0 | blue == 1),
         green_zero_one = (green == 0 | green == 1))

df.ind = data.pe %>%
  dplyr::select(blue, green, AC, `A-C`, `-AC`, `-A-C`, prolific_id, id) %>% 
  mutate(bg_ind = blue * green,
         diff = AC - bg_ind,
         min_bg = case_when(blue + green > 1 ~ blue + green - 1, T ~ 0),
         max_bg = pmin(blue, green), 
         min_is_max_bg = round(min_bg, 2) == round(max_bg, 2),
         diff_ac_min = round(AC, 2) - round(min_bg, 2),
         diff_ac_max = round(max_bg, 2) - round(AC, 2)
         ) %>% 
  filter(str_detect(id, "independent"))

ind_trials = df.ind$id %>% unique()

# normalize rating if necessary, that is, when AC is (due to rounding) not 
# in necessary possible range, e.g. blue: 0.8, green: 0.8 --> min: 0.6, but
# observed value for AC is 0.59
df.to_normalize = df.ind %>% 
  mutate(ac_out_range = AC < min_bg | AC > max_bg) %>% 
  filter(ac_out_range & !min_is_max_bg & diff_ac_min !=0 & diff_ac_max != 0)
  # check how often observed is limit of possible interval:
  # filter(diff_ac_max == 0 | diff_ac_min == 0)
  # not in range and observed ac is not same as min and max value for ac:

df.normalized = df.to_normalize %>% 
  group_by(id, prolific_id) %>% 
  mutate(s = sum(AC, `A-C`, `-AC`, `-A-C`), 
         AC = AC/s,
         `A-C` = `A-C`/s, 
         `-AC` = `-AC`/s, 
         `-A-C` = `-A-C`/s, 
         blue = AC + `A-C`, 
         green = AC + `-AC`, 
         bg_ind = blue * green,
         diff = AC - bg_ind,
         min_bg = case_when(blue + green > 1 ~ blue + green - 1, T ~ 0),
         max_bg = pmin(blue, green)
  )

# all data with those normalized where necessary
df.behav <- bind_rows(
  anti_join(df.ind, df.normalized %>% dplyr::select(prolific_id, id)),
  df.normalized %>%
    dplyr::select(-s, -ac_out_range)
) %>% 
  mutate(
    max_sub = bg_ind - min_bg, # at most this subtracted from perfect indep P(a,c)=P(b)*P(g)
    max_add = max_bg - bg_ind # at most this added to perfect indep P(a,c)=P(b)*P(g)
  )

# Plot data ---------------------------------------------------------------
zoib_plots = function(df.behav){
  df.behav.long <- df.behav %>% pivot_longer(cols=c(blue, green)) 
  df.behav.beta <- df.behav.long %>% filter(!value %in% c(0, 1))
  p.beta = df.behav.beta %>%
    ggplot(aes(x=value)) + geom_density(aes(color=name)) +
    facet_wrap(~id, scales="free") +
    scale_color_manual(values = c("blue"="darkblue", "green"="forestgreen"))
  
  df.behav.zoi <- df.behav.long  %>% filter(value %in% c(0,1)) %>% 
    mutate(value = as.factor(value)) %>% 
    group_by(id, name, value) %>%
    summarize(n = n(), .groups = "drop_last") %>%
    mutate(proportion = n / sum(n))
  
  p.zoi = df.behav.zoi %>%
    ggplot(aes(x=name, y = proportion)) +
    geom_bar(stat="identity", aes(fill=value), position=position_dodge()) +
    facet_wrap(~id) + labs(y = "P(x=value|zoi)", x = "probability")
  return(list(zoi=p.zoi, beta=p.beta))
}
zoib_plots(df.behav)


# Some tests --------------------------------------------------------------
# check ratio difference P(b,g)-P(b)*P(g)
df.behav %>%
  ggplot(aes(x = diff)) +
  geom_density() +
  facet_wrap(~id, scales = "free")
# Simulate some data ------------------------------------------------------
N=89
fake_data = tibble(
  id = "independent_hl", 
  blue = simulate_zoib_data(N, alpha=0.3, gamma=0.92, 
                           shape1=6, shape2=2),
  green = simulate_zoib_data(N, alpha=0.4, gamma=0.98, 
                            shape1=2, shape2=6),
  delta = rnorm(N, mean=0, sd=0.05)
  ) %>% 
  mutate(bg_ind = blue * green,
         min_bg = case_when(blue + green > 1 ~ blue + green - 1, T ~ 0),
         max_bg = pmin(blue, green),
         max_sub = bg_ind - min_bg, # at most this subtracted from perfect indep P(a,c)=P(b)*P(g)
         max_add = max_bg - bg_ind, # at most this added to perfect indep P(a,c)=P(b)*P(g)
         diff = case_when(delta < 0 & delta < -1 * max_sub ~ max_sub,
                          delta < 0 ~ delta,
                          delta > 0 & delta > max_sub ~ max_sub,
                          delta > 0 ~ -delta),
         AC = bg_ind + diff,
         ) %>% 
  add_column(prolific_id = "")
# plot simulated fake data  
zoib_plots(fake_data)
fake_data %>% ggplot(aes(x=delta)) + geom_density()
fake_data %>% ggplot(aes(x=diff)) + geom_density()

# Run webppl --------------------------------------------------------------
posterior_samples.ind = map_dfr(ind_trials, function(trial_id){
  message(trial_id)
  df.trial = df.behav %>% filter(id == trial_id)
  data_webppl = list(probs = df.trial)
  samples.posterior <-
    posterior_samples.ind <- webppl(
    program_file = here("webppl-model", "posterior-independent-data.wppl"),
    data_var = "data",
    model_var = "non_normalized_posterior",
    data = data_webppl,
    packages = c(paste("webppl-model", "webppl-packages", "dataHelpers", sep = FS)),
    inference_opts = list(method = "incrementalMH",
                          samples = 5000,
                          lag = 10,
                          burn = 50000,
                          verbose = F),
    chains = 4, 
    cores = 4) %>%
    as_tibble() %>% 
    pivot_wider(names_from = "Parameter", values_from = "value") %>% 
    add_column(id = trial_id)
  
  return(samples.posterior)
})
save_data(posterior_samples.ind, paste(target_dir, "posterior_samples.rds", sep=FS))
# posterior_samples.ind <- readRDS(paste(target_dir, "posterior_samples.rds", sep=FS))
posterior_samples.long <- posterior_samples.ind %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_longer(cols=c(-Iteration, -Chain, -id), names_to = "param")

# Chain Plots
plots_chains = map(c("alpha", "gamma", "shape1", "shape2"), function(par){
  plots = map(c("blue", "green"), function(p) {
    fn <- paste(par, p, sep="_")
    p <- posterior_samples.long %>% 
      filter(startsWith(param, fn)) %>% 
      ggplot(aes(x=Iteration, y = value, color=Chain, group = Chain)) + 
      geom_line() +
      theme(axis.text = element_text(size = 18)) +
      facet_wrap(id~param, scales="free", 
                 labeller = labeller(id = trial_names))
    ggsave(paste(target_dir, FS, paste("chains_", fn, ".png", sep = ""), sep=""), 
           p, width = 15, height = 11)
    
    p2 <- posterior_samples.long %>% 
      filter(startsWith(param, fn)) %>% 
      ggplot(aes(x=value, color=Chain, group = Chain)) + 
      geom_density() + 
      theme(axis.text = element_text(size = 18)) +
      facet_wrap(id~param, scales="free", 
                 labeller = labeller(id = trial_names))
    ggsave(paste(target_dir, FS, 
                 paste("chains_marginal_", fn, ".png", sep = ""), sep=""), p2,
           width = 15, height = 11)
    return(p)                 
  })
  return(plots)
})
# compute posterior means and MCMC diagnostics
df.diagnostics = map(ind_trials, function(trial_id){
  zoib_data = map(c("alpha", "gamma", "shape1", "shape2"), function(par){
    vals = map(c("blue", "green"), function(p) {
      fn <- paste(par, p, sep="_")
      samples <- posterior_samples.long %>% 
        filter(param == !!fn & id == !!trial_id)
      mat <- samples %>% dplyr::select(Iteration, Chain, value) %>% 
        rename(`.iteration`=Iteration, `.chain` = Chain) %>% 
        posterior::as_draws_matrix()
      df.summary = posterior::summarise_draws(mat) %>% mutate(variable = fn)
      return(df.summary %>% add_column(id = trial_id))
    })
  }) %>% bind_rows()
  # sigma
  samples <- posterior_samples.long %>% filter(param=="sd_delta" & id==trial_id)
  mat <- samples %>% dplyr::select(Iteration, Chain, value) %>% 
    rename(`.iteration`=Iteration, `.chain` = Chain) %>% 
    posterior::as_draws_matrix()
  df.summary = posterior::summarise_draws(mat) %>%
    mutate(variable = "sd_delta", id=trial_id)
  return(bind_rows(df.summary, zoib_data))
  }) %>% bind_rows()

save_data(df.diagnostics, paste(target_dir, "mcmc-diagnostics.rds", sep=FS))
#df.diagnostics <- readRDS(paste(target_dir, "mcmc-diagnostics.rds", sep=FS))

df.diagnostics %>% filter(rhat >=1.1)


# Mean posterior values for likelihood function
posterior_means <- df.diagnostics %>% dplyr::select(variable, mean, id) %>% 
  group_by(id) %>% 
  separate(variable, into=c("par", "prob"), sep="_") %>% 
  arrange(prob)
posterior_means %>% filter(id == "independent_uh")
write_csv(posterior_means %>% 
            mutate(mean = round(mean, 2)) %>% 
            pivot_wider(names_from = "par", values_from = "mean"), 
          paste(target_dir, FS, "posterior_means_likelihood_fns.csv", sep=""))

# Plots with expected values of posterior distributions
plots.pp_evs = map(ind_trials, function(trial_id){
  df.trial = df.behav %>% filter(id == !!trial_id)
  pars = posterior_means %>% filter(id == !!trial_id) %>% 
    unite("name", par, prob, sep="_") %>% 
    pivot_wider(names_from = "name", values_from="mean")
  N <- nrow(df.trial)
  mat = map(seq(1, 100), function(i){
    sim.blue = simulate_zoib_data(N, pars$alpha_blue, 
                                  pars$gamma_blue,
                                  pars$shape1_blue, 
                                  pars$shape2_blue)
    sim.green = simulate_zoib_data(N, pars$alpha_green, 
                                   pars$gamma_green,
                                   pars$shape1_green, 
                                   pars$shape2_green)
    sim.diff = rnorm(N, mean = 0, sd = pars$sd_delta)
    return(c(y=c(sim.blue, sim.green, sim.diff)))
  }) %>% bind_rows() %>% as.matrix()
  grp <- factor(c(rep("blue", N), rep("green", N), rep("delta", N)), 
                levels = c("blue", "green", "delta"))
  y <- bind_rows(c(y=c(df.trial$blue, df.trial$green, df.trial$diff)))
  
  p <- ppc_dens_overlay_grouped(y = y[1,] %>% as.numeric(), 
                                yrep = mat, group = grp) +
    labs(title = trial_names[[trial_id]]) +
    theme(axis.text.x = element_text(size = 18), 
          panel.spacing = unit(2, "lines")) +
    facet_wrap("group", scales = "free")
  fn <- paste(target_dir, paste("new_data_evs_", trial_id, ".png", sep=""), sep=FS)
  ggsave(fn, p, width=11)
  return(p)
})
plots.pp_evs


# Generate new independent tables -----------------------------------------
# generate set of independent tables for each sample from posterior distribution (parameters)
sampled_tables <- sample_tables(df.behav, ind_trials, 
                                posterior_samples.ind %>% filter(Iteration <= 25),
  here("webppl-model", "posterior-independent-data.wppl")
)
save_data(sampled_tables, 
          paste(target_dir, "sampled-tables-posterior-predictive.rds", sep=FS))
# sampled_tables <- readRDS(paste(target_dir,
#                                 "sampled-tables-posterior-predictive.rds", 
#                                 sep=FS))

# sampled_tables.evs <- sample_tables(
#   df.behav, ind_trials, 
#   posterior_means %>% 
#     unite("name", par, prob, sep="_") %>% 
#     pivot_wider(names_from = "name", values_from = "mean"),
#   here("webppl-model", "posterior-independent-data.wppl"),
#   repetitions = 20
# )

# plot new tables sampled for each independent context with observed data
make_pp_plots = function(df.behav, trials, sampled_tables, target_dir, fn_prefix){
  pp_plots <- map(trials, function(trial_id){
    df.trial <- df.behav %>% filter(id == trial_id)
    df.samples <- sampled_tables %>% filter(id == trial_id)
    tit <- trial_names[[trial_id]]
    
    fn <- paste(target_dir, FS, paste(fn_prefix, "_", trial_id, ".png", sep=""), sep="")
    p <- plot_new_tables(df.trial, df.samples, tit, fn)
  })
  return(pp_plots)
}
pp_plots <- make_pp_plots(df.behav, ind_trials, sampled_tables, target_dir, "pp-tables")


# Posterior predictive ----------------------------------------------------
# Log likelihood plots
likelihood_fn = function(df.samples, df.grp){
  message(paste(df.grp$id))
  data_webppl <- list(probs = df.behav %>% filter(id == df.grp$id),
                      samples_posterior = df.samples)
  samples.pp <- webppl(
    program_file = here("webppl-model", "posterior-predictive-independent-data.wppl"),
    random_seed = params$seed_webppl,
    data_var = "data",
    data = data_webppl,
    packages = c(paste("webppl-model", "webppl-packages", "dataHelpers", sep = FS))
  ) 
  
  result = tibble(
    x_blue = samples.pp$x_blue,
    x_green = samples.pp$x_green,
    x_bg = samples.pp$x_bg,
    ll_X_new = samples.pp$ll_X, 
    ll_obs = samples.pp$ll_obs, 
    ll_X_blue = samples.pp$ll_x_blue, 
    ll_X_green = samples.pp$ll_x_green, 
    ll_X_bg = samples.pp$ll_x_bg, 
    ll_obs_blue = samples.pp$ll_obs_blue, 
    ll_obs_green = samples.pp$ll_obs_green, 
    ll_obs_bg = samples.pp$ll_obs_bg
  )
  return(result %>% add_column(id = df.grp$id))
}
pp_samples_ll = group_map(posterior_samples.ind %>% group_by(id),
                          likelihood_fn) %>% bind_rows()
save_data(pp_samples_ll, paste(target_dir, "pp_samples_ll.rds",sep = FS))
# pp_samples_ll <- readRDS(paste(target_dir, "pp_samples_ll.rds",sep = FS))

# overall log likelihood of data from posterior predictive
ll_X_obs.mean = pp_samples_ll %>% group_by(id) %>% 
  summarize(ev = mean(ll_obs), .groups = "keep")

p <- pp_samples_ll %>% 
  ggplot(aes(x = ll_X_new)) + geom_density() +
  facet_wrap(~id, scales = "free", labeller = labeller(id = trial_names)) +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  labs(x = "log likelihood", y = "density")
p
ggsave(paste(target_dir, "pp-log-likelihood-independent.png", sep = FS), p)

# log likelihood separate for P(b), P(g) and P(b,g)
pp_samples_ll.long <- pp_samples_ll %>%
  pivot_longer(cols = c(-id,-x_blue, -x_green, -x_bg), 
               names_to = "prob", names_prefix = "ll_", values_to = "ll")
# ll pp-samples
df.ll_X = pp_samples_ll.long %>% 
  filter(startsWith(prob, "X_") & prob != "X_new") %>% 
  mutate(prob = str_replace(prob, "X_", ""))

# ev ll observed data given params from posterior
ll_X_obs.mean = pp_samples_ll.long %>% group_by(id, prob) %>% 
  filter(startsWith(prob, "obs_")) %>% 
  summarize(ev = mean(ll), .groups = "drop_last") %>% 
  mutate(prob = str_replace(prob, "obs_", ""))

p <- df.ll_X %>% 
  ggplot(aes(x = ll)) + geom_density() +
  facet_wrap(prob~id, ncol=5, scales = "free", 
             labeller = labeller(id = trial_names)) +
  geom_point(data = ll_X_obs.mean, aes(x=ev, y=0), size=2, color = 'firebrick') +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18)
        ) + 
  labs(x = "log likelihood", y = "density")
p
ggsave(paste(target_dir, "pp-log-likelihood-independent-separate.png", sep=FS), 
       p, height = 9, width = 18)

