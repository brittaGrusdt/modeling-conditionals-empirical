library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
# Priors
active_config = "default_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
if(!dir.exists(params$dir_results)) dir.create(params$dir_results, recursive = T)
result_dir = params$dir_results

# behavioral data
data.behav <- read_csv(here(params$dir_data, params$fn_cleaned_data)) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task, slider) %>% 
  translate_standardized2model() 

# add parameters to run webppl model
# observed data each participant and trial + ratios
pars.observations <- get_observed_data(data.behav, params, cutoff=F)

# draw samples from prior
pars.rsa_states <- get_rsa_states(params)
# load likelihood parameters fitted to data
pars.likelihoods <- get_likelihood_params_fitted_data(params)
# if not provided, all utterances equally likely
# params$p_utts = rep(1 / length(params$utterances), length(params$utterances))

weights <- readRDS(
  here(params$dir_results,
       paste("weights_ci_", params$n_forward_samples, ".rds", sep=""))
) %>% filter(prior_r == "informative")

p_s_ci = list(p_s_ci = left_join(weights, pars.rsa_states$prior_samples, by="bn_id"))

Sys.setenv(R_CONFIG_ACTIVE = "priors_relations")
params$prior_relations <- config::get()[["informative"]]

params <- c(params, pars.observations, pars.likelihoods, pars.rsa_states, p_s_ci)
params$par_fit <- c("alpha", "theta", "utt_cost") #, "gamma")
params$packages <- c(params$packages, 
                     paste("webppl-model", "node_modules", "dataHelpers", sep = FS))
# set to utterance cost
# params$utt_cost <- params_evs %>% 
#   filter(!Parameter %in% c("alpha", "theta")) %>% 
#   pivot_wider(names_from = "Parameter", values_from = "value")

mcmc_params <- tibble(n_samples = 10000, n_burn = 0, n_lag = 5, n_chains = 4)
posterior <- webppl(program_file = here("webppl-model", "fit-rsa.wppl"), 
                    data_var = "data",
                    model_var = "non_normalized_posterior",
                    data = params,
                    inference_opts = list(method = "MCMC",
                                          samples = mcmc_params$n_samples,
                                          burn = mcmc_params$n_burn,
                                          lag = mcmc_params$n_lag,
                                          verbose = T),
                    packages = params$packages, 
                    chains = mcmc_params$n_chains, 
                    cores = mcmc_params$n_chains
                    ) %>% as_tibble()

posterior_samples <- posterior %>% unnest(c(value)) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  add_column(mcmc = list(mcmc_params), n_forward_samples = params$n_forward_samples)

fn <- str_flatten(params$par_fit, collapse = "_")
save_data(posterior, 
          here(params$dir_results, 
               paste("mcmc-posterior-fit-context-predictions-", fn, ".rds", sep="")
               ))
               
# chain plot, iteration vs. value for each chain
p_chain = posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  ggplot(aes(x=Iteration, y = value, color = Chain)) +
  geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_both, ncol = 3) 
p_chain
ggsave(here(params$dir_results, paste("chain_", fn, ".png", sep="")), p_chain)


# utterance prior probabilities (cost)
posterior_utts <- posterior_samples %>%
  filter(startsWith(as.character(Parameter), "utts.")) %>% 
  mutate(Parameter = as.character(Parameter), 
         Parameter = str_replace(Parameter, "utts.", ""), 
         utterance = Parameter) %>% 
  chunk_utterances() %>% rename(utt_type = utterance) %>% 
  mutate(Parameter = factor(Parameter, 
                            levels = c("C and A", "C and -A", "-C and A", "-C and -A",
                                       "A", "C", "-A", "-C",
                                       "might A", "might C", "might -A", "might -C",
                                       "A > C", "A > -C", "C > A", "C > -A", 
                                       "-A > C", "-A > -C", "-C > A", "-C > -A" 
                                       )))
posterior_utts %>% 
  ggplot(aes(x=value, color = Chain)) +
  geom_density() + 
  facet_wrap(~Parameter, scales = "free_y", labeller = label_value, ncol = 4) 


# expected values
params_evs <- posterior_samples %>% group_by(Parameter) %>% 
  filter(Iteration >= 1000 & Chain == 3) %>% 
  summarize(value = mean(value), .groups = "drop_last") %>%
  mutate(Parameter=str_replace(Parameter, "utts.", ""))
params_evs

# softmax utt cost
softmax <- function(vec,i){return(exp(vec[i])/sum(exp(vec)))}
costs <- params_evs %>% filter(!Parameter %in% c("alpha", "theta", "gamma")) %>%
  mutate(value = -1 * value) %>% pull(value)
df.p_utts = params_evs %>% filter(!Parameter %in% c("alpha", "theta", "gamma")) %>% 
  add_column(p = softmax(costs)) %>% arrange(desc(p))

params_evs %>% pivot_wider(names_from = "Parameter", values_from = "value") 



posterior_utts %>% group_by(Chain, Parameter) %>% 
  summarize(mean = mean(value)) %>%
  arrange(desc(mean)) %>% 
  pivot_wider(names_from = "Parameter", values_from = "mean")

posterior_utts %>% 
  ggplot(aes(x=Iteration, y = value, color = Chain)) +
  geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_value, ncol = 5) 


# plot posterior densities
posterior_samples %>% filter(Parameter %in% c("alpha", "theta")) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  ggplot(aes(x=value, color = Chain)) +
  geom_density() + 
  facet_wrap(Parameter~Chain, scales = "free", labeller = label_both, ncol = 4) 

# alpha vs. theta
posterior_samples %>% filter(Parameter %in% c("alpha", "theta")) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_wider(names_from="Parameter", values_from = "value") %>% 
  ggplot(aes(x=alpha, y=theta, color = Chain)) +
  geom_point()








