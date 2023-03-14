library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))

theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)

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

# add parameters to run webppl model
# observed data each participant and trial + ratios
pars.observations <- get_observed_data(data.behav, params)
# draw samples from prior
pars.rsa_states <- get_rsa_states(params)
# sanity check:
pars.rsa_states$prior_samples %>% 
  group_by(probability) %>% dplyr::count() %>% ungroup() %>%
  mutate(proportion = n / sum(n), N=sum(n))

# load likelihood parameters fitted to data
pars.likelihoods <- get_likelihood_params_fitted_data(params)
# if not provided, all utterances equally likely
# params$p_utts = rep(1 / length(params$utterances), length(params$utterances))
params <- c(params, pars.observations, pars.likelihoods, pars.rsa_states)

# Run Model ---------------------------------------------------------------
# 1. predictions by contexts
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()
params$prior_relations <- par_relations[["uninformative"]]
params$utt_cost <- tibble(utterance = params$utterances, cost=params$utt_cost)

params$alpha <- 1.98
params$theta <- 0.764

params$alpha <- 2.58
params$theta <- 0.687
# params$gamma <- 1
# params$utt_cost <- df.p_utts %>% dplyr::select(-p) %>% 
#   pivot_wider(names_from="Parameter", values_from="value")

model <- "var rsa_predictions = run_rsa_model(data)
rsa_predictions
"
data <-   webppl(program_code = model,
                 data = params,
                 data_var = "data",
                 random_seed = params$seed_webppl,
                 packages = params$packages
)
posterior <- data %>% map(function(x){as_tibble(x)})

###############################################################################
# new prediction based on repeatedly drawn N_participants states from 
# weights to make rsa predictions
model.predictions <- posterior %>% bind_rows() %>% rowid_to_column() %>% 
  unnest(c(p_hat, utterance)) %>% 
  mutate(ll_ci = as.numeric(ll_ci)) %>% 
  arrange(id, p_hat)
#hdis.mean = mean_hdi(model.predictions %>% group_by(id, utterance), p_hat)
pred.quantiles = group_map(model.predictions %>% group_by(id, utterance), 
                           function(df, df.grp){
  qs = quantile(df$p_hat, c(0.025, 0.5, 0.975))
  tibble(estimate=qs[['50%']], lower=qs[['2.5%']], upper = qs[['97.5%']], 
         id=df.grp$id, 
         utterance=df.grp$utterance)
}) %>% bind_rows()
model.predictions %>% filter(str_detect(id, "if1_uh")) %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1") +
  labs(x = "predicted probability")
###############################################################################

model.predictions <- posterior %>% 
  map(function(x){as_tibble(x) %>% mutate(ll_ci = as.numeric(ll_ci))}) %>% 
  bind_rows() %>% group_by(id) %>% 
  mutate(p_hat_round = round(p_hat, 2)) %>% 
  arrange(desc(p_hat))


production.joint = left_join(model.predictions, 
                             params$observed_utts_ratios %>% rename(p = ratio)) %>% 
  rename(model = p_hat, behavioral = p) %>% 
  mutate(relation = case_when(startsWith(id, "independent") ~ "independent", 
                              startsWith(id, "if1") ~ "if1",
                              startsWith(id, "if2") ~ "if2"))

plot_correlation(production.joint, color = "id") + facet_wrap(~id)
plot_model_vs_data_bars(production.joint, 
                        str_flatten(c("alpha", params$alpha, "theta", 
                                      params$theta), collapse="_"), 
                        by_utt_type = F)













