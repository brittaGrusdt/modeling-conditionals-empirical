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
source(here("R", "helpers-rsa-model.R"))
theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
config_cns = "fine_grained_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "semi_informative"
config_fits <- "rsa_fit_params"
config_speaker_type <- "pragmatic_utt_type"
params <- prepare_data_for_wppl(config_cns, config_weights_relations, 
                                config_fits = config_fits,
                                config_speaker_type = config_speaker_type,
                                extra_packages = extra_packages)

# empirically observed data
df.bootstrapped_uc_ratios_ci <- readRDS(params$bootstrapped_ci_ratios) %>% 
  translate_standardized2model() %>% 
  mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                  utts.model.literals,
                                                  utts.model.ifs,
                                                  utts.model.mights)))
obs.ind <- df.bootstrapped_uc_ratios_ci %>%
  filter(str_detect(id, "independent")) %>% 
  mutate(trial=id, id = map_chr(id, get_name_context)) %>% 
  rename(utt = utterance) %>% 
  mutate(utterance = utt.standardized) %>% 
  chunk_utterances() %>% 
  rename(utt_type = utterance, utterance = utt)
obs.dep <- df.bootstrapped_uc_ratios_ci %>%
  filter(str_detect(id, "if")) %>% mutate(trial = id) %>% 
  rename(utt = utterance) %>% 
  mutate(utterance = utt.standardized) %>% 
  chunk_utterances() %>% 
  rename(utt_type = utterance, utterance = utt)
  
# Prior predictive --------------------------------------------------------
# condition statements are necessary to make sure that we can actually run 
# the RSA-model with the sampled set of parameters!

model <- "
var model_prior = function(){
  var params = priorSample(data['par_fit'])
  setParams(params)
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
          paste(params$config_dir, "samples-prior.rds", sep=FS))

# chain plot
prior_samples %>% filter(!startsWith(Parameter, "utts.")) %>% 
  ggplot(aes(x=Iteration, y=value, color=Chain)) + geom_line() +
  facet_wrap(~Parameter, scales="free")

# density plot with MCMC-samples from prior distributions overlayed with 
# samples from underlying prior distributions (without considering invalid
# utt-states based on theta)
x <- seq(0, 100, by=0.01)
distrs = tibble(y_alpha = exp(rnorm(length(x), mean=params$alpha_mu, sd=params$alpha_sigma)),
                y_theta = rbeta(length(x), params$theta_shape1, params$theta_might_shape2)) %>% 
  pivot_longer(cols = c(y_alpha, y_theta), names_to = "Parameter", names_prefix = "y_")
# check different parametrizations for alpha # 
distrs %>% ggplot(aes(x=value)) + geom_density() + xlim(0, 300)
  
prior_samples %>% 
  filter(!startsWith(Parameter, "utts.")) %>% 
  ggplot(aes(x=value)) + 
  geom_density(aes(color=Chain)) +
  geom_density(data=distrs) +
  facet_wrap(~Parameter, scales="free") 

prior_samples %>% pivot_wider(names_from="Parameter", values_from="value") %>% 
  ggplot(aes(x=alpha, y=theta)) + 
  geom_density_2d_filled()

# then run RSA-model once with each sampled set of parameters
params$sampled_params <- format_param_samples(prior_samples)
params$verbose <- 0
rsa_data <- webppl(program_file = params$wppl_predictive_checks,
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
save_data(prior_predictive, 
          paste(params$config_dir, "prior-predictive.rds", sep=FS))
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
prior_pred.ind %>% #filter(id=="ind:HH") %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1") +
  labs(x = "predicted probability")

# 2. prior predictive bar plot for each utterance dependent contexts if1
prior_pred.dep %>% filter(str_detect(id, "if1")) %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
  facet_wrap(~utterance, scales="free_y", ncol=4) +
  theme(axis.text.x = element_text(size=8)) +
  scale_fill_brewer(name = "context", palette = "Set1", 
                    labels = labels_dep_contexts) +
  labs(x = "predicted probability") 

# 3. prior predictive bar plot for each utterance dependent contexts if2
prior_pred.dep %>% filter(str_detect(id, "if2")) %>% 
  ggplot(aes(x=p_hat, fill=id)) + 
  geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
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


# Prior Predictive in dependence of alpha and theta -----------------------
ranges_alpha <- c("0 < alpha < 1",  "1 <= alpha < 2",
                  "2 <= alpha < 5", "5 <= alpha < 10",
                  "10 <= alpha < 20", "20 <= alpha < 30", "30 <= alpha")
pp <- left_join(
  prior_predictive, 
  params$sampled_params %>% dplyr::select(alpha, theta, sample_id)
)  %>% mutate(utt=utterance, utterance = as.character(utterance)) %>% 
  chunk_utterances() %>% 
  rename(utt_type = utterance, utterance = utt) %>% 
  mutate(cat_alpha = case_when(alpha < 1 ~ ranges_alpha[1], 
                               alpha < 2 ~ ranges_alpha[2], 
                               alpha < 5 ~ ranges_alpha[3], 
                               alpha < 10 ~ ranges_alpha[4],
                               alpha < 20 ~ ranges_alpha[5],
                               alpha < 30 ~ ranges_alpha[6],
                               T ~ ranges_alpha[7]),
         cat_alpha = factor(cat_alpha, levels = ranges_alpha)
  )
pp %>% group_by(cat_alpha) %>% dplyr::count() %>% ungroup() %>% 
  mutate(ratio = n/sum(n))

trial <- "independent_hh"
trial <- "if2_ll"
ut <- "conjunction"
if(str_detect(trial, "independent")){
  df.obs <- obs.ind %>% filter(trial == !!trial & utt_type == ut)
}else {
  df.obs <- obs.dep %>% filter(trial == !!trial & utt_type == ut)
}
pp %>% filter(id==trial & utt_type == ut) %>%  
  ggplot(aes(color = utterance, group=utterance)) + 
  geom_point(aes(x=theta, y = p_hat), alpha=0.5) +
  geom_line(aes(x=theta, y = p_hat), alpha=0.5) +
  geom_point(data = df.obs, aes(y=estimate, x=0), size=1, stroke=1, shape=4) +
  # geom_segment(data = df.obs, aes(x=0, y=lower, xend = 0.8, yend=lower)) +
  # geom_segment(data = df.obs, aes(x=0, y=upper, xend = 0.8, yend=upper)) +
  facet_wrap(~cat_alpha) +
  ggtitle(get_name_context(trial))



# Log likelihoods ---------------------------------------------------------
pp.ll = pp %>% distinct_at(vars(c(id, alpha, theta, ll_ci))) 
# only for independent conrtexts ll becomes -Inf
pp.ll %>% filter(is.infinite(ll_ci)) %>% pull(id) %>% unique()

pp.ll %>% filter(str_detect(id, "independent") & !is.infinite(ll_ci)) %>% 
  arrange(desc(theta))
min_theta_inf = pp.ll %>% 
  filter(str_detect(id, "independent") & is.infinite(ll_ci)) %>% 
  group_by(id) %>% filter(theta == min(theta)) 
alpha = min_theta_inf$alpha %>% unique
theta = min_theta_inf$theta %>% unique

pp %>% 
  filter(alpha== !!alpha & str_detect(id, "independent")) %>% 
  filter(theta >= !!theta) %>% 
  arrange(id, p_hat) %>% 
  filter(p_hat == 0)

# are there any dependent contexts where any utterance is predicted with probab. 0?
# NO
pp %>% dplyr::select(id, p_hat, utterance) %>% filter(str_detect(id, "if")) %>% 
  filter(p_hat == 0)

# check weights for independent contexts
weights.assertable_u = params$weights %>% get_assertable_utterances(theta)

# weights for problematic utterance 'neither block falls' are very small 
# for dependent contexts (given large theta so that -Inf in indep. contexts)
# and =0 for independent contexts!
weights.assertable_u %>% group_by(id) %>% 
  filter(utt == standardized.sentences$none) %>% 
  distinct_at(vars(c(id, bn_id)), .keep_all = T) %>% 
  summarize(summed_weight = sum(probs)) %>% arrange(desc(summed_weight))


# log likelihood plot
pp.ll %>%
  ggplot(aes(x=theta, y=ll_ci, color=id, group=id)) +
  geom_point() + geom_line() + labs(y="log likelihood data") +
  facet_wrap(~id, scales="free")


# Posterior predictive ----------------------------------------------------
subfolder <- create_subconfig_folder_for_fitting(
  config_dir = params$config_dir, 
  par_fit = params$par_fit, 
  speaker_type = params$speaker_type
)
# get samples from posterior distribution
posterior_samples <- readRDS(here(params$config_dir, subfolder, "mcmc-posterior.rds"))
# just run rsa once for each identical parameter combination!
params$sampled_params <- format_param_samples(posterior_samples) %>% 
  distinct_at(vars(c(alpha, theta)), .keep_all = T)

evs <- posterior_samples %>% group_by(Parameter) %>% 
  summarize(value = mean(value))
# then run RSA-model once with each sampled set of parameters
data <- webppl(program_file = params$wppl_predictive_checks,
               data = params,
               data_var = "data",
               random_seed = params$seed_webppl,
               packages = params$packages)


# TODO: add predictions as often as param combi occurred
posterior_predictive <- data %>% imap(function(x, id){
  predictions <- as_tibble(x)
  df.predictions <- map(predictions, function(y){
    as_tibble(y) %>% mutate(ll_ci = as.numeric(ll_ci))
  }) %>% bind_rows() %>% add_column(sample_id = id)
}) %>% bind_rows()
save_data(posterior_predictive, 
          paste(params$config_dir, subfolder, "posterior-predictive.rds", sep=FS))

predictives <- bind_rows(
  prior_predictive %>% add_column(distribution = "prior predictive"),
  posterior_predictive %>% add_column(distribution = "posterior predictive")
)

trials <- params$observations$id %>% unique()
utterances <- posterior_predictive$utterance %>% unique()
map(trials, function(trial){
  if(str_detect(trial, "independent")){
    df.obs <- obs.ind %>% filter(trial == !!trial)
  }else {
      df.obs <- obs.dep %>% filter(trial == !!trial)
  }
  # add 0 for unobserved utterances
  df.obs <- bind_rows(
    df.obs, 
    tibble(utterance = utterances[!utterances %in% df.obs$utterance],
          estimate = 0, lower = 0, upper = 0, trial = trial, 
          id = df.obs$id[1])
  ) %>% 
  mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                  utts.model.literals,
                                                  utts.model.ifs,
                                                  utts.model.mights)))
  
  p <- predictives %>%
    mutate(utterance = factor(utterance, levels = c(utts.model.conjs, 
                                                    utts.model.literals,
                                                    utts.model.ifs,
                                                    utts.model.mights))) %>% 
    filter(str_detect(id, trial)) %>% 
    ggplot() +
    geom_histogram(aes(x=p_hat, fill=distribution), alpha=0.5,
                   position = 'identity', binwidth=0.005) +
    geom_point(data = df.obs, aes(x=estimate, y=0), color='firebrick', 
               size=1, stroke=1, shape=4) +
    geom_segment(data = df.obs, aes(x=lower, y=0, xend = upper, yend=0), 
                 color = 'firebrick') +
    #geom_line(data = df.obs, aes(x=lower, y=upper, group = utterance), size=1) +
    facet_wrap(~utterance, scales="free", ncol=4) +
    theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8)) +
    scale_fill_brewer(name = "distribution", palette = "Set1") +
    labs(x = "predicted probability", y="count", title = get_name_context(trial))
  ggsave(paste(params$config_dir, FS, subfolder, FS, 
               "predictive-checks-", trial, ".png", sep=""), 
         p, width = 10, height = 10)
})


# Predictions for single states -------------------------------------------
# to see how alpha and theta influence predictions for a particular state, 
# e.g. where a conjunction is assertable

# params$utt_cost <- tibble(utterance = params$utterances, cost=params$utt_cost)

wppl_code = "
  var sampled_params = data['sampled_params']
  var bn_ids = _.map(ALL_BNS, 'bn_id')
  var idx_bn = bn_ids.indexOf(data['state_id'][0])
  var bn = ALL_BNS[idx_bn]
  var s_id = bn.bn_id
  // iterate over sampled parameters
  var rsa_predictions = map(function(pars){
    var params = {
      alpha: pars.alpha,
      //gamma: pars.gamma,
      theta: pars.theta,
      //utt_cost: Object.fromEntries(
      //zip(_.map(data.utt_cost, 'utterance'), _.map(data.utt_cost, 'cost'))
      //)
    }
    display('theta: ' + params.theta + ' alpha: ' + params.alpha)
    setParams(params)
  
    var rsa_speaker = run_speaker([bn], data['speaker_type'], false)
    var state_prob_pairs = Object.values(Object.values(rsa_speaker[s_id].params)[0])
    var utts = _.map(state_prob_pairs, 'val')
    var ps = _.map(state_prob_pairs, 'prob')
    
    return({p_hat: ps,
            utterance: utts, 
            theta: params.theta,
            alpha: params.alpha,
            sample_id: pars.sample_id,
            bn_id: s_id})
  }, sampled_params)

  rsa_predictions
"

# single state where P(b,g) is closest to 1
states <- params$prior_samples %>%
  dplyr::select(bn_id, r, probability, table.probs, table.support) %>%
  unnest(c(table.probs, table.support)) %>%
  group_by(bn_id) %>%
  pivot_wider(names_from = "table.support", values_from = "table.probs")
s = states %>% ungroup() %>% arrange(desc(AC)) %>% slice(0, 1)

theta=0.799
utts.assertable <- states %>% get_assertable_utterances(theta) %>% 
  filter(bn_id == s$bn_id) %>% pull(utt)

model_pars = format_param_samples(prior_samples) %>% arrange(theta)
params$sampled_params <- model_pars[seq(1, nrow(model_pars), by=10),]
params$state_id = s$bn_id

#keep alpha at 1
#params$sampled_params$alpha <- 1
predictions.speaker = webppl(
  program_code = wppl_code,
  data = params,
  data_var = "data",
  packages = params$packages
) %>% as_tibble() %>% unnest(c(p_hat, utterance)) %>% 
  mutate(utt=utterance) %>%  chunk_utterances() %>% 
  rename(utt_type = utterance, utterance = utt) %>% 
  mutate(cat_alpha = case_when(alpha < 1 ~ "0 < alpha < 1", 
                               alpha < 2 ~ "1 <= alpha < 2", 
                               alpha < 5 ~ "2 <= alpha < 5", 
                               alpha < 10 ~ "5 <= alpha < 10",
                               T ~ "10 <= alpha", 
  ))

predictions.speaker %>% group_by(cat_alpha) %>% dplyr::count() %>% ungroup() %>% 
  mutate(ratio = n/sum(n))

predictions.speaker %>% 
  ggplot(aes(x=theta, y = p_hat, color = utt_type, group=utt_type)) + 
  geom_point() + geom_line() +
  facet_wrap(~cat_alpha)


