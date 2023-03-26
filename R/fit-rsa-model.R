library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)
library(ggdist)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))
theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
config_cns = "fine_grained_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "flat_dependent"
config_fits <- "alpha_theta_gamma"
config_speaker_type <- "pragmatic_utt_type"

params <- prepare_data_for_wppl(config_cns, config_weights_relations,
                                config_fits = config_fits,
                                config_speaker_type = config_speaker_type,
                                extra_packages = extra_packages)
mcmc_params <- tibble(n_samples = 2500, n_burn = 2500, n_lag = 5, n_chains = 4)
posterior <- webppl(program_file = params$wppl_fit_rsa, 
                    data_var = "data",
                    model_var = "non_normalized_posterior",
                    data = params,
                    inference_opts = list(method = "incrementalMH",
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
  add_column(mcmc = list(mcmc_params), nb_rsa_states = params$nb_rsa_states)

save_data(posterior_samples %>% 
            add_column(config_prior_r = config_weights_relations, 
                       config_cns = config_cns), 
          here(params$speaker_subfolder, "mcmc-posterior.rds"))
# posterior_samples <- readRDS(here(params$speaker_subfolder, "mcmc-posterior.rds")) 
              
# chain plot, iteration vs. value for each chain
p_chain = posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  ggplot(aes(x=Iteration, y = value, color = Chain)) +
  geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_both, ncol = 3) 
p_chain
ggsave(here(params$speaker_subfolder, "chains.png"), p_chain)

# plot posterior densities
p.density_posterior = posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  ggplot(aes(x=value, color = Chain)) +
  geom_density() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_parsed, ncol = 4)
p.density_posterior
ggsave(here(params$speaker_subfolder, "density_posterior.png"), p.density_posterior)

# posterior with highest density intervals
hdis.mean_posterior <- mean_hdi(posterior_samples %>% group_by(Parameter), value)
hdis.mean_posterior

p.posterior_hdis = posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  ggplot(aes(x = value, fill = Parameter)) +
  stat_halfeye(.width = c(0.95), point_interval = "median_hdi") +
  facet_wrap(~Parameter, labeller = label_parsed, scales = "free") +
  theme(legend.position = "none") +
  labs(x="posterior value", y = "density")
p.posterior_hdis
ggsave(here(params$speaker_subfolder, "density_posterior_hdis.png"), p.posterior_hdis)

# alpha vs. theta
posterior_samples %>% filter(Parameter %in% c("alpha", "theta")) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_wider(names_from="Parameter", values_from = "value") %>% 
  ggplot(aes(x=alpha, y=theta, color = Chain)) +
  geom_point()
# gamma vs. theta
posterior_samples %>% filter(Parameter %in% c("gamma", "theta")) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  pivot_wider(names_from="Parameter", values_from = "value") %>% 
  ggplot(aes(x=gamma, y=theta, color = Chain)) +
  geom_point(alpha=0.5)

# expected values
params_evs <- posterior_samples %>% group_by(Parameter) %>% 
  summarize(value = mean(value), .groups = "drop_last") %>%
  mutate(Parameter=str_replace(Parameter, "utts.", ""))
params_evs %>% pivot_wider(names_from = "Parameter", values_from = "value") 

# plots for fitted utterance cost -----------------------------------------
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

# softmax utt cost
softmax <- function(vec,i){return(exp(vec[i])/sum(exp(vec)))}
costs <- params_evs %>% filter(!Parameter %in% c("alpha", "theta", "gamma")) %>%
  mutate(value = -1 * value) %>% pull(value)
df.p_utts = params_evs %>% filter(!Parameter %in% c("alpha", "theta", "gamma")) %>% 
  add_column(p = softmax(costs)) %>% arrange(desc(p))

posterior_utts %>% group_by(Chain, Parameter) %>% 
  summarize(mean = mean(value)) %>%
  arrange(desc(mean)) %>% 
  pivot_wider(names_from = "Parameter", values_from = "mean")

posterior_utts %>% 
  ggplot(aes(x=Iteration, y = value, color = Chain)) +
  geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_value, ncol = 5) 





