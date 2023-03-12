library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)
library(tidybayes)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

# parameters to run webppl model
params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
result_dir = params$dir_results

path_model_file = paste(params$dir_wppl_code, "predictive-checks.wppl", sep=FS)

################################################################################
# use more fine-grained causal nets (causal power, noise: low/unc/high)
active_config = "fine_grained_causal_nets"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
pars.cns <- config::get()
cns = create_causal_nets(pars.cns$dep_causal_power, pars.cns$dep_noise, 
                         pars.cns$dep_marginal, params$rels_dep,
                         pars.cns$ind_a, pars.cns$ind_c)
params$causal_nets_dep = cns$dep
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
# load likelihood parameters fitted to data
pars.likelihoods <- get_likelihood_params_fitted_data(params)
# if not provided, all utterances equally likely
# params$p_utts = rep(1 / length(params$utterances), length(params$utterances))
params <- c(params, pars.observations, pars.likelihoods, pars.rsa_states)

# empirically observed data
df.observed <- params$observed_utts_ratios %>% 
  mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                  utts.model.literals,
                                                  utts.model.ifs,
                                                  utts.model.mights)))
obs.ind <- df.observed %>% filter(str_detect(id, "independent")) %>% 
  mutate(id = map_chr(id, get_name_context))
obs.dep <- df.observed %>% filter(str_detect(id, "if")) 
  
# Run Model ---------------------------------------------------------------
# 1. predictions by contexts
params$packages <- c(params$packages, paste("webppl-model", "node_modules", 
                                            "dataHelpers", sep = FS))

active_config = "priors_relations"
Sys.setenv(R_CONFIG_ACTIVE = active_config)
par_relations <- config::get()
params$prior_relations <- par_relations[["informative"]]
#params$utt_cost <- tibble(utterance = params$utterances, cost=params$utt_cost)

# Helper functions --------------------------------------------------------
format_param_samples = function(samples){
  samples.utt_cost <- samples %>% 
    filter(startsWith(as.character(Parameter), "utts.")) %>% 
    mutate(Parameter = as.character(Parameter), 
           Parameter = str_replace(Parameter, "utts.", "")) %>% 
    rename(cost = value, utterance = Parameter) %>% 
    group_by(Iteration, Chain) %>% 
    nest() %>% rename(utt_cost=data)
  
  samples.formatted <- left_join(
    samples.utt_cost,
    samples %>% 
      filter(!startsWith(as.character(Parameter), "utts.")) %>% 
      pivot_wider(names_from=Parameter, values_from=value)
  ) %>% 
    mutate(sample_id = paste("chain", Chain, "-Iteration", Iteration, sep=""))
  return(samples.formatted)
}

# Prior predictive --------------------------------------------------------
params$par_fit <- c("alpha", "theta", "utt_cost") #, "gamma")

# condition statements are necessary to make sure that we can actually run 
# the RSA-model with the sampled set of parameters!
model <- "
var model_prior = function(){
  var params = priorSample(data['par_fit'])
  setParams(params)
  //display(params)
  // with current sampled set of parameters, are there any states where no
  // utterance is applicable or any utterance without state?
  // if yes do not consider these parameters
  var invalid_utts_states = check_states_utts(
    ALL_BNS,
    _.map(globalStore.utterances, 'utt'),
    globalStore.thresholds,
    params,
    false
  )
  condition(invalid_utts_states.states.length == 0)
  condition(invalid_utts_states.utts.length == 0)

  query.add('utts', params.utt_cost)
  query.add('alpha', params.alpha)
  query.add('theta', params.theta)
  query.add('gamma', params.gamma)
  return(query)
}
"

mcmc_params <- tibble(n_samples = 1000, n_burn = 1000, n_lag = 15, n_chains = 4)
# get samples from prior distribution
prior_samples <- webppl(
  program_code = model,
  data = params,
  data_var = "data",
  model_var = "model_prior",
  inference_opts = list(method = "incrementalMH",
                        samples = mcmc_params$n_samples,
                        burn = mcmc_params$n_burn,
                        lag = mcmc_params$n_lag,
                        verbose = T),
  packages = params$packages, 
  chains = mcmc_params$n_chains, 
  cores = mcmc_params$n_chains
  ) %>% as_tibble() %>% 
  mutate(Parameter = as.character(Parameter), 
         Parameter = str_replace(Parameter, "utt_cost", "utts"), 
         Chain = as.factor(Chain))
save_data(prior_samples %>% add_column(n_burn = mcmc_params$n_burn, 
                                       n_lag = mcmc_params$n_lag), 
          paste(result_dir, "samples-prior.rds", sep=FS))

# chain plot
prior_samples %>% filter(!startsWith(Parameter, "utts.")) %>% 
  ggplot(aes(x=Iteration, y=value, color=Chain)) + geom_line() +
  facet_wrap(~Parameter, scales="free")

# density plot with MCMC-samples from prior distributions overlayed with 
# samples from underlying prior distributions (without considering invalid
# utt-states based on theta)
x <- seq(0, 200, by=0.01)
distrs = tibble(y_alpha = exp(rnorm(length(x), mean=1.5, sd=1)),
                y_theta = rbeta(length(x), 8, 2)) %>% 
  pivot_longer(cols = c(y_alpha, y_theta), names_to = "Parameter", names_prefix = "y_")

prior_samples %>% 
  filter(!startsWith(Parameter, "utts.")) %>% 
  ggplot(aes(x=value)) + 
  geom_density(aes(color=Chain)) +
  geom_density(data=distrs) +
  facet_wrap(~Parameter, scales="free") 
# distribution of theta is shifted away from 1 since for larger theta there are
# utterances without states!


# then run RSA-model once with each sampled set of parameters
params$sampled_params <- format_param_samples(prior_samples)[1:100,]
rsa_data <- webppl(program_file = path_model_file,
                   data = params,
                   data_var = "data",
                   random_seed = params$seed_webppl,
                   packages = params$packages)

prior_predictive <- rsa_data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows() %>% 
  mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                  utts.model.literals,
                                                  utts.model.ifs,
                                                  utts.model.mights)))
# Plots -------------------------------------------------------------------
model_utts = c("-A", "A", "-C", "C",
               "-C and -A", "-C and A", "C and -A", "C and A", 
               "might -A", "might A", "might -C", "might C", 
               "A > C", "C > A", "A > -C", "-C > A",
               "-A > C", "C > -A", "-A > -C", "-C > -A"
)
model_utt_labels =  c("¬b", "b", "¬g", "g",
                      "¬b¬g", "b¬g", "¬bg", "bg", 
                      "might ¬b", "might b", "might ¬g", "might g", 
                      "if b,g", "if g,b", "if b,¬g", "if ¬g,b",
                      "if ¬b,g", "if g, ¬b", "if ¬b,¬g", "if ¬g, ¬b"
)

medians.prior <- median_hdi(prior_predictive %>% group_by(id, utterance), p_hat) %>% 
  distinct_at(vars(c(id, utterance, p_hat)))
means.prior <- mean_hdi(prior_predictive %>% group_by(id, utterance), p_hat) %>%
  distinct_at(vars(c(id, utterance, p_hat)))

prior_pred.ind <- prior_predictive %>% filter(str_detect(id, "independent")) %>% 
  mutate(id = map_chr(id, get_name_context))
prior_pred.dep <- prior_predictive %>% filter(str_detect(id, "if"))

# 1. prior predictive bar plot for each utterance independent contexts
prior_pred.ind %>% #filter(id=="ind:UH") %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  # geom_point(data = obs.ind, aes(x=ratio, y=0, color=id), size=1, stroke=1, shape=4) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1") +
  labs(x = "predicted probability")

# 2. prior predictive bar plot for each utterance dependent contexts if1
prior_pred.dep %>% filter(str_detect(id, "if1")) %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  # geom_point(data = obs.dep %>% filter(str_detect(id, "if1")), 
  #            aes(x=ratio, y=0, color=id), size=1, stroke=1, shape=4) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1", 
                    labels = labels_dep_contexts) +
  labs(x = "predicted probability") 

# 3. prior predictive bar plot for each utterance dependent contexts if2
prior_pred.dep %>% filter(str_detect(id, "if2")) %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  # geom_point(data = obs.dep %>% filter(str_detect(id, "if2")), 
  #            aes(x=ratio, y=0, color=id), size=1, stroke=1, shape=4) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1", 
                    labels = labels_dep_contexts) +
  labs(x = "predicted probability")


# Density plot for utterances for a single context
trial <- "if1_uh"
prior_predictive %>% filter(id==trial) %>%
  ggplot(aes(x=p_hat, fill=utterance, color=utterance)) +
  geom_density() +
  stat_halfeye(color='darkgrey', .width = c(0.95), point_interval = "mean_hdi") +
  theme(legend.position = "none", axis.text.x = element_text(size=10)) +
  facet_wrap(~utterance, scales="free", ncol = 4) +
  ggtitle(get_name_context(trial))

# Predictions based on mean/median of parameters (from prior distributions)
means.prior %>% mutate(p_hat = round(p_hat, 3)) %>% filter(p_hat > 0) %>% 
  mutate(utt=utterance, utterance=as.character(utterance)) %>%
  chunk_utterances() %>% 
  rename(utt_type=utterance, utterance=utt) %>% 
  mutate(utterance = factor(utterance, levels = model_utts, labels = model_utt_labels)) %>% 
  filter(str_detect(id, "independent")) %>% 
  ggplot(aes(y=utterance, x = p_hat, fill=utt_type)) + 
  geom_bar(stat="identity") +
  theme(axis.text.y = element_text(size=10)) +
  facet_wrap(~id, scales = "free")



# Posterior predictive ----------------------------------------------------
fn <- "alpha_theta_utt_cost"
# get samples from posterior distribution
posterior_samples <- readRDS(here(params$dir_results, 
                             paste("mcmc-posterior-fit-context-predictions-",
                                   fn, ".rds", sep="")))
params$sampled_params <- format_param_samples(posterior_samples)[1:100,]

# then run RSA-model once with each sampled set of parameters
data <- webppl(program_file = path_model_file,
               data = params,
               data_var = "data",
               random_seed = params$seed_webppl,
               packages = params$packages)

posterior_predictive <- data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows()


predictives <- bind_rows(
  prior_predictive %>% add_column(distribution = "prior predictive"),
  posterior_predictive %>% add_column(distribution = "posterior predictive")
)
predictives %>%
  mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                  utts.model.literals,
                                                  utts.model.ifs,
                                                  utts.model.mights))) %>% 
  filter(str_detect(id, trial)) %>% 
  ggplot(aes(x=p_hat)) + 
  geom_histogram(alpha=0.75, aes(fill=distribution)) +
  geom_point(data = df.observed, aes(x=ratio, y=0), size=2, color='black') +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8))



