library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(ggthemes)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))

theme_set(theme_clean(base_size = 26) + 
            theme(legend.position = "top", text = element_text(size = 26)))
# Setup -------------------------------------------------------------------
config_cns = "fine_grained_dep_cns"
extra_packages = c("dataHelpers")
config_weights_relations = "flat_dependent"
config_speaker_type = "pragmatic_utt_type"
# config_speaker_type = "literal"
config_fits = "alpha_theta_gamma"
# config_fits = "alpha_theta"
params <- prepare_data_for_wppl(config_cns = config_cns, 
                                config_weights_relations = config_weights_relations, 
                                config_speaker_type = config_speaker_type,
                                config_fits = config_fits,
                                extra_packages = extra_packages)
speaker_subfolder <- params$speaker_subfolder
# set alpha and theta, otherwise default values used
# params$alpha <- 2.29
# params$theta <- 0.859
# params$gamma <- 0.296
# params$utt_cost <- df.p_utts %>% dplyr::select(-p) %>% 
#   pivot_wider(names_from="Parameter", values_from="value")

# use MAP-parameters
speaker_model <- str_split(config_speaker_type, "_")[[1]][1]
if(str_ends(config_fits, "_gamma")) {
  speaker_model <- paste(speaker_model, "gamma", sep="_")
}

params <- add_MAP_params(params, speaker_model, config_speaker_type)

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

# save predictions using MAP-parameters
save_data(model.predictions %>% add_column(theta = params$theta,
                                           gamma = params$gamma,
                                           alpha = params$alpha),
          here(params$speaker_subfolder, paste("rsa-results-MAP.rds", sep=FS)))

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

p.corr_ind = plot_correlation(
  production.joint %>% filter(relation == "independent"), 
  ncol = 3
)
p.corr_ind
ggsave(paste(params$speaker_subfolder, "corr-plot-single-run_IND_MAPs.png", 
             sep = FS),
       p.corr_ind, width = 26, height = 12)

p.corr_dep = plot_correlation(
  production.joint %>% filter(relation != "independent"), 
  ncol = 4
)
p.corr_dep
ggsave(paste(params$speaker_subfolder, "corr-plot-single-run_DEP_MAPS.png", 
             sep=FS), 
       p.corr_dep, width = 30, height = 15)

# plot_model_vs_data_bars(production.joint, 
#                         str_flatten(c("alpha", params$alpha, "theta", 
#                                       params$theta), collapse="_"), 
#                         by_utt_type = F)
