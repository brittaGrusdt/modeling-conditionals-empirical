library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)
library(ggdist)
library(ggthemes)
library(bayesplot)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))
theme_set(theme_clean(base_size = 26) + 
            theme(legend.position = "top", text = element_text(size = 26)))

# Setup -------------------------------------------------------------------
config_cns = "fine_grained_dep_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "flat_dependent"
config_speaker_type <- "pragmatic_utt_type"
#config_speaker_type <- "literal"
config_fits <- "alpha_theta"
params <- prepare_data_for_wppl(config_cns, config_weights_relations,
                                config_fits = config_fits,
                                config_speaker_type = config_speaker_type,
                                extra_packages = extra_packages, 
                                mcmc_params = "mcmc_default_params")
posterior <- webppl(program_file = params$wppl_fit_rsa, 
                    data_var = "data",
                    model_var = "non_normalized_posterior",
                    data = params,
                    inference_opts = list(method = "incrementalMH",
                                          samples = params$mcmc$n_samples,
                                          burn = params$mcmc$n_burn,
                                          lag = params$mcmc$n_lag,
                                          verbose = T),
                    packages = params$packages, 
                    chains = params$mcmc$n_chains, 
                    cores = params$mcmc$n_chains
                    ) %>% as_tibble()

posterior_samples <- posterior %>% unnest(c(value)) %>% 
  mutate(Chain = as.factor(Chain)) %>% 
  add_column(mcmc = list(params$mcmc), nb_rsa_states = params$nb_rsa_states)

save_data(posterior_samples %>% 
            add_column(config_prior_r = config_weights_relations, 
                       config_cns = config_cns), 
          paste(params$speaker_subfolder, "mcmc-posterior.rds", sep = FS))
# posterior_samples <- readRDS(here(params$speaker_subfolder, "mcmc-posterior.rds")) 


if(str_detect(config_fits, "gamma")){
  pars <- c("theta", "gamma")
} else {
  pars <- c("theta")
}
if(config_speaker_type != "literal"){
  pars <- c(pars, "alpha")
} 
# check pairs plot
df.posterior <- posterior_samples %>% dplyr::select(Chain, Iteration, Parameter, value) %>% 
  pivot_wider(names_from = "Parameter", values_from = "value", values_fn = list) %>% 
  unnest(all_of(pars)) %>%
  mutate(Chain = as.integer(Chain))

if(str_detect(config_fits, "gamma") || config_speaker_type != "literal"){
  color_scheme_set("blue")
  p.pairs <- mcmc_pairs(df.posterior %>% dplyr::select(-Iteration))
  ggsave(here(params$speaker_subfolder, "pairs.png"), p.pairs)
}
  
df.diagnostics <- df.posterior %>%  
  rename(`.iteration`=Iteration, `.chain` = Chain) %>% 
  posterior::as_draws_matrix() %>% 
  posterior::summarise_draws()
df.diagnostics
write_csv(df.diagnostics, here(params$speaker_subfolder, "diagnostics.csv"))

# chain plot, iteration vs. value for each chain
p.chain <- mcmc_trace(df.posterior %>% dplyr::select(-Iteration), 
                      facet_args = list(labeller = as_labeller(
                        x = c('alpha' = 'alpha', 'theta' = 'theta', 
                              'gamma' = 'gamma'), 
                        default = label_parsed)
                        )) + theme(legend.position = "top") 
p.chain
p_chain = posterior_samples %>% 
  filter(!startsWith(as.character(Parameter), "utts.")) %>% 
  ggplot(aes(x=Iteration, y = value, color = Chain)) +
  geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_parsed, ncol = 3) 
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

p.posterior_hdis <- 
  mcmc_dens(df.posterior, pars = pars, 
            facet_args = list(labeller = as_labeller(
              x = c('alpha' = 'alpha', 'theta' = 'theta', 'gamma' = 'gamma'),
              default = label_parsed
              ))) +
  geom_segment(data = hdis.mean_posterior,
               aes(x=.lower, xend = .upper, y=0, yend=0), 
               color = 'red', size = 2) +
  geom_point(data=hdis.mean_posterior, aes(x=value, y=0), color='black', size=2.5) +
  theme(panel.spacing = unit(2, "lines"))
p.posterior_hdis  
ggsave(here(params$speaker_subfolder, "density_posterior_hdis.png"), p.posterior_hdis)

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





