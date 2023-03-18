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
config_weights_relations = "semi_informative"
params <- prepare_data_for_wppl(config_cns, config_weights_relations, 
                                extra_packages = extra_packages)

# set alpha and theta, otherwise default values used
params$alpha <- 3.37
params$theta <- 0.338
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

###############################################################################
# new prediction based on repeatedly drawn N_participants states from 
# weights to make rsa predictions
# model.predictions <- posterior %>% bind_rows() %>% rowid_to_column() %>% 
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













