library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))

theme_set(theme_minimal(base_size=20) + theme(legend.position = "top"))

# Setup -------------------------------------------------------------------
config_cns = "fine_grained_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "flat_dependent"
config_speaker_type = "pragmatic_utt_type"
params <- prepare_data_for_wppl(config_cns = config_cns, 
                                config_weights_relations = config_weights_relations, 
                                config_speaker_type = config_speaker_type,
                                extra_packages = extra_packages)
# set alpha and theta, otherwise default values used
# params$alpha <- 2.29
# params$theta <- 0.859
# params$gamma <- 0.296
# params$utt_cost <- df.p_utts %>% dplyr::select(-p) %>% 
#   pivot_wider(names_from="Parameter", values_from="value")

# use MAP-parameters
speaker_model <- "pragmatic_gamma"
if(str_split(speaker_model, "_")[[1]][1] != 
   str_split(config_speaker_type, "_")[[1]][1]) stop("mismatch speaker model and speaker type")

Sys.setenv(R_CONFIG_ACTIVE = paste("MAP", speaker_model, sep="_"))
pars.speaker_model <- config::get()
params$alpha <- pars.speaker_model$alpha
params$theta <- pars.speaker_model$theta
params$gamma <- pars.speaker_model$gamma


model <- "var rsa_predictions = run_rsa_model(data)
rsa_predictions
"
data <-   webppl(program_code = model,
                 data = params,
                 data_var = "data",
                 random_seed = params$seed_webppl,
                 packages = params$packages
)
wppl_output <- data %>% map(function(x){as_tibble(x)})
model.predictions <- wppl_output %>% 
  map(function(x){as_tibble(x) %>% mutate(ll_ci = as.numeric(ll_ci))}) %>% 
  bind_rows() %>% group_by(id) %>% 
  mutate(p_hat_round = round(p_hat, 2)) %>% 
  arrange(desc(p_hat))

# save_data(model.predictions %>% add_column(theta=params$theta,
#                                            gamma = params$gamma,
#                                            alpha = params$alpha),
#           here(params$config_dir,
#                paste("rsa-results-MAP-", speaker_model, ".rds", sep="")))

model.predictions %>% dplyr::select(ll_ci, id) %>% distinct() %>% 
  arrange(desc(ll_ci))

production.joint = left_join(model.predictions, 
                             params$observed_utts_ratios %>% rename(p = ratio)) %>% 
  rename(model = p_hat, behavioral = p) %>% 
  mutate(relation = case_when(startsWith(id, "independent") ~ "independent", 
                              startsWith(id, "if1") ~ "if1",
                              startsWith(id, "if2") ~ "if2")) %>% 
  mutate(label_id = map_chr(id, get_str_contexts)) %>% 
  mutate(response = utterance) %>% 
  translate_utterances() %>% 
  rename(utt = utterance, utterance = response)

p.corr = plot_correlation(production.joint)
ggsave(paste(params$config_dir, FS,
             "alpha_theta_gamma", FS,
             params$speaker_type, FS, 
             "corr-plot-single-run_alpha-", params$alpha, 
             "-theta-", params$theta, 
             "-gamma-", params$gamma, 
             ".png", sep=""), 
       p.corr, 
       width = 21, height = 10)
plot_model_vs_data_bars(production.joint, 
                        str_flatten(c("alpha", params$alpha, "theta", 
                                      params$theta), collapse="_"), 
                        by_utt_type = F)

###############################################################################
# new prediction based on repeatedly drawn N_participants states from 
# weights to make rsa predictions
# model.predictions <-  wppl_output %>% bind_rows() %>% rowid_to_column() %>% 
#   unnest(c(p_hat, utterance)) %>% 
#   mutate(ll_ci = as.numeric(ll_ci)) %>% 
#   arrange(id, p_hat)
# #hdis.mean = mean_hdi(model.predictions %>% group_by(id, utterance), p_hat)
# pred.quantiles = group_map(model.predictions %>% group_by(id, utterance), 
#                            function(df, df.grp){
#   qs = quantile(df$p_hat, c(0.025, 0.5, 0.975))
#   tibble(estimate=qs[['50%']], lower=qs[['2.5%']], upper = qs[['97.5%']], 
#          id=df.grp$id, 
#          utterance=df.grp$utterance)
# }) %>% bind_rows()
# model.predictions %>% filter(str_detect(id, "if1_uh")) %>% 
#   ggplot(aes(x=p_hat, fill=id)) + 
#   geom_histogram(alpha=0.5, position = 'identity', binwidth=0.01) +
#   facet_wrap(~utterance, scales="free_y", ncol=4) +
#   theme(axis.text.x = element_text(size=8)) +
#   scale_fill_brewer(name = "context", palette = "Set1") +
#   labs(x = "predicted probability")
###############################################################################













