library(here)
library(rwebppl)
library(tidyverse)
library(tibble)
library(ExpDataWrangling)
library(ModelUtils)
library(stringr)
library(tidybayes)
library(ggpubr)
library(ggdist)
library(ggthemes)
library(scales)

source(here("R", "helpers-plotting.R"))
source(here("R", "helpers-load-data.R"))
source(here("R", "helpers-rsa-model.R"))
theme_set(theme_clean(base_size = 26) + 
            theme(legend.position = "top", text = element_text(size = 26)))

# Setup
config_weights_relations <- "flat_dependent"
config_cns = "fine_grained_dep_cns"
extra_packages = c("dataHelpers")
# use pragmatic extended speaker model here for calling fn to load all data later
config_speaker_type = "pragmatic_utt_type"
config_fits <- "alpha_theta_gamma"
params <- prepare_data_for_wppl(config_cns = config_cns,
                                config_weights_relations = config_weights_relations, 
                                extra_packages = extra_packages,
                                config_speaker_type =  config_speaker_type,
                                config_fits = config_fits
)
subfolder_pragmatic_sp_extended <- params$speaker_subfolder



# empirical data
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

df.observed <- bind_rows(obs.ind, obs.dep)
TRIALS <- params$observations$id %>% unique()

# Plots across speaker types ----------------------------------------------
# Log-likelihood plots 
map.ll.speakers <- get_data_all_speakers(
  config_weights_relations, config_cns, extra_packages, 
  "MAP-log-likelihoods.csv", 
  params$speaker_subfolder
) %>% group_by(alpha, theta, gamma, speaker_model)

# plot MAP ll for each context with overall MAP-parameters
# summed ll across contexts
map.ll.speakers.summed <- map.ll.speakers %>%
  group_by(sample_id, speaker_model, alpha, theta, gamma) %>% 
  summarize(summed_ll = round(sum(ll_ci)), .groups = "drop_last") %>%  
  #summed_neg_ll = round(sum(neg_ll_ci)), .) %>% 
  mutate(alpha = round(alpha, 2), 
         theta = round(theta, 2),
         gamma = case_when(!is.na(gamma) ~ round(gamma, 2), 
                           T ~ gamma), 
         label_par = case_when(speaker_model == "literal" ~ paste("theta:", theta), 
                               speaker_model == "literal.gamma" ~ 
                                 paste("theta: ", theta, " gamma: ", gamma, sep=""),
                               speaker_model == "pragmatic" ~ 
                                 paste("alpha: ", alpha, " theta: ", theta, sep=""),
                               speaker_model == "pragmatic.gamma" ~ 
                                 paste("alpha: ", alpha, " theta: ", theta, 
                                     " gamma: ", gamma, sep=""),
                               speaker_model == "random" ~ ""),
         y = case_when(speaker_model == "random" ~ -175,
                       speaker_model == "literal.gamma" ~ -155, 
                       speaker_model == "literal" ~ -165, 
                       speaker_model == "pragmatic" ~ -90, 
                       speaker_model == "pragmatic.gamma" ~ -40)
  )
map.ll.speakers.summed

p.MAP_ll <- map.ll.speakers %>% 
  mutate(label_id = map_chr(id, get_str_contexts)) %>% 
  ggplot(aes(x=label_id, y = ll_ci, color = speaker_model, group = speaker_model)) + 
  geom_point() + geom_line() +
  geom_text(data = map.ll.speakers.summed %>% add_column(x="ind:UL"), 
            aes(x=x, y=y, label = summed_ll), size = 6) +
  geom_text(data = map.ll.speakers.summed %>% add_column(x="'if'[2]*':HL'"), 
            aes(x=x, y=y, label = label_par), size = 6) +
  labs(x = "context", y = "Log likelihood\n for MAP parameters") +
  scale_x_discrete(labels = label_parse()) +
  scale_color_manual(name = "speaker model", 
                     values = SPEAKER_COLORS[names(SPEAKER_COLORS) != "observed"])
p.MAP_ll
ggsave(filename = here(params$mcmc_folder, "MAP_ll.png"), p.MAP_ll)

# Boxplots log likelihoods for each sample from posterior (alpha, theta, gamma)
pp.ll.speakers <- get_data_all_speakers(
  config_weights_relations, 
  config_cns, 
  extra_packages,
  "posterior_predictive_ll.rds", 
  params$speaker_subfolder
)
p.pp_ll_speakers <- pp.ll.speakers %>% 
  # filter(!endsWith(speaker_model, ".gamma")) %>% 
  mutate(label_id = map_chr(id, get_str_contexts)) %>% 
  ggplot(aes(x=label_id, y=ll_ci, color=speaker_model)) + 
  geom_boxplot() +
  scale_x_discrete(labels = label_parse()) +
  scale_color_manual(name = "speaker model", 
                     values = SPEAKER_COLORS[names(SPEAKER_COLORS) != "observed"]) +
  labs(x="context", y = "log likelihood")
p.pp_ll_speakers
ggsave(filename = here(params$mcmc_folder, "pp_ll_speakers.png"), 
       p.pp_ll_speakers, width = 17, height = 8)

# MAP-values for each context
MAP.contexts <- get_data_all_speakers(
  config_weights_relations, 
  config_cns, 
  extra_packages, 
  "MAP-log-likelihood-ci.csv", 
  params$speaker_subfolder
)
p.MAPs_contexts = MAP.contexts %>% filter(speaker_model != "random") %>% 
  mutate(label_id = map_chr(id, get_str_contexts)) %>% 
  pivot_longer(cols = c(alpha, theta, gamma), names_to = "Parameter") %>% 
  filter(!is.na(value)) %>% # gamma is NA if non-existent
  # mutate(value = case_when(is.na(value) ~ 1, T ~ value)) %>% 
  ggplot(aes(x = label_id, y = value, color = speaker_model, group = speaker_model)) + 
  geom_point() + geom_line() + 
  facet_wrap(~Parameter, scales = "free", labeller = label_parsed, ncol = 2) + 
  theme(axis.text.x = element_text(size=11)) +
  scale_x_discrete(labels = label_parse()) + 
  scale_color_manual(name = "speaker model", 
                     values = SPEAKER_COLORS[c("literal", "pragmatic", 
                                               "literal.gamma", "pragmatic.gamma")]) +
  theme(axis.text.x = element_text(size = 16)) + 
  labs(x = "context", y = "MAP parameter value")
p.MAPs_contexts
ggsave(filename = here(params$mcmc_folder, "MAPs_contexts.png"), p.MAPs_contexts, 
       height = 10, width = 22)


# Data vs. Model with MAP-parameters --------------------------------------
# run RSA model with MAP-parameters (across all contexts) for each speaker model
# (using script 'single-model-run.R')

path_pragmatic_ext <- paste(subfolder_pragmatic_sp_extended, 
                            "rsa-results-MAP.rds", sep=FS)
path_pragmatic <- str_replace(path_pragmatic_ext, "alpha_theta_gamma", 
                              "alpha_theta")
predictions.MAP <- bind_rows(
  readRDS(path_pragmatic_ext) %>% add_column(speaker_model = "pragmatic.gamma"),
  readRDS(str_replace(path_pragmatic_ext, "pragmatic_utt_type", "literal")) %>% 
    add_column(speaker_model = "literal.gamma"),
  readRDS(path_pragmatic) %>% add_column(speaker_model = "pragmatic"),
  readRDS(str_replace(path_pragmatic, "pragmatic_utt_type", "literal")) %>% 
    add_column(speaker_model = "literal")
)

joint <- left_join(
  predictions.MAP,
  df.observed %>% dplyr::select(estimate, lower, upper, utterance, trial) %>% 
    rename(id = trial)
) %>% 
  mutate(estimate = case_when(is.na(estimate) ~ 0, T ~ estimate),
         lower = case_when(is.na(lower) ~ 0, T ~ lower),
         upper = case_when(is.na(upper) ~ 0, T ~ upper)) %>% 
  mutate(response = utterance) %>% 
  translate_utterances() %>% 
  mutate(label_par = 
           case_when(is.na(alpha) ~ "",
                     gamma == 1 ~ paste("alpha: ", alpha, " theta: ", theta, sep=""),
                     T ~ paste("alpha: ", alpha, " theta: ", theta, " gamma: ", gamma, sep="")))

plots <- group_map(joint %>% group_by(id), function(df, df.grp){
  pars <- df %>% dplyr::select(label_par, speaker_model, alpha, theta, gamma) %>% distinct() 
  p <- df %>% 
    ggscatter(y = "estimate", x = "p_hat", add = "reg.line",
              conf.int = TRUE, cor.method = "pearson",
              ylab = "Data", xlab = "Prediction (MAP parameters)") +
    # geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
    geom_point(aes(x=p_hat, y=estimate, color = response, fill = response)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color=response)) +
    geom_text(data = pars, 
              aes(x=Inf, y=Inf, label = paste("alpha:", alpha, 
                                              "~''*gamma:", gamma,
                                              "~''*theta:", theta)),
              size = 6, parse=T, vjust = 1.5, hjust = 1) +
    facet_wrap(~speaker_model, scales = "free") +
    scale_color_manual(name = "utterance", values = UTT_COLORS) +
    scale_fill_manual(name = "utterance", values = UTT_COLORS) +
    guides(fill = guide_legend(title.position = "top"), 
           color = guide_legend(title.position = "top")) +
    theme(text = element_text(size=16), legend.text = element_text(size=14)) +
    stat_cor(size = 6) +
    ggtitle(get_name_context(df.grp$id))
  ggsave(paste(params$mcmc_folder, FS, "models-vs-data-MAP-", df.grp$id, ".png", 
               sep=""), p, width = 16, height = 9)
  return(p)
})

# check neg log likelihood ind:UL -----------------------------------------
levels_utts <- c("-A", "A", "-C", "C", 
                 "-C and -A", "-C and A", "C and -A",  "C and A",
                 "might -A",  "might A",   "might -C",  "might C", 
                 "A > C", "C > A", "A > -C", "-C > A", 
                 "-A > C", "C > -A", "-A > -C",   "-C > -A")
pred.ind_ul <- predictions.MAP %>% 
  filter(speaker_model == "literal" & id=="independent_ul") %>% 
  mutate(utterance = factor(utterance, levels = levels_utts)) %>% 
  arrange(utterance)

pred.ind_ul$p_hat


# other plots -------------------------------------------------------------
# group relevant conditionals together
context <- "if2_hl"
joint.grp <- joint %>% 
  mutate(utt.grp = case_when(utterance %in% c("A > C", "C > A", "-A > -C", "-C > -A") ~ "conditional", 
                             T ~ response)) %>% 
  filter(id==!!context & speaker_model == "pragmatic") %>% 
  group_by(utt.grp) %>% 
  summarize(estimate = sum(estimate), p_hat = sum(p_hat))

pars <- joint %>% 
  filter(id==!!context & speaker_model == "pragmatic") %>% 
  dplyr::select(alpha, theta, gamma, label_par) %>% distinct() 

# adapt utterance colors
utt_colors <- c(`both blocks fall` = "darkgoldenrod1", 
                `blue falls but green does not fall` = "brown1", 
                `green falls but blue does not fall` = "brown3",
                `neither block falls` = "darkred", 
                `blue falls` = "skyblue1",
                `green falls` = "green",
                `blue does not fall` = "blue", 
                `green does not fall` = "darkgreen",
                `conditional` = "orchid1",
                `if blue falls green does not fall` = "gray15",
                `if green falls blue does not fall` = "gray30",
                `if blue does not fall green falls` = "gray50",
                `if green does not fall blue falls` = "gray80", 
                `blue might fall` = "cyan",
                `green might fall` = "yellow2",
                `blue might not fall` = "turquoise4",
                `green might not fall` = "yellow4")


p <- joint.grp %>% 
  ggscatter(y = "estimate", x = "p_hat", add = "reg.line",
            conf.int = TRUE, cor.method = "pearson",
            ylab = "Data", xlab = "Prediction (MAP parameters)") +
  geom_point(aes(x=p_hat, y=estimate, color = utt.grp, fill = utt.grp)) +
  geom_text(data = pars, 
            aes(x=Inf, y=Inf, label = paste("alpha:", alpha, 
                                            "~''*gamma:", gamma,
                                            "~''*theta:", theta)),
            size = 6, parse=T, vjust = 1.5, hjust = 1) +
  guides(fill = guide_legend(title.position = "top"), 
         color = guide_legend(title.position = "top")) +
  theme(text = element_text(size=16), legend.text = element_text(size=14)) +
  scale_color_manual(name = "utterance", values = utt_colors) +
  scale_fill_manual(name = "utterance", values = utt_colors) +
  stat_cor(size = 6) +
  ggtitle(get_name_context(context))
p


# Some more checks --------------------------------------------------------
# use of conditionals
df.use_conditionals <- params$observed_utts_ratios %>% 
  filter(str_detect(utterance, ">") & n > 0) %>% 
  arrange(desc(n)) %>% 
  mutate(conditional.bg = startsWith(utterance, "A") | startsWith(utterance, "-A"),
         conditional.gb = startsWith(utterance, "C") | startsWith(utterance, "-C")
  ) 

df.use_conditionals %>% 
  group_by(conditional.bg, conditional.gb) %>% 
  summarize(n = sum(n))


df.use_conditionals %>% filter(conditional.gb)
df.use_conditionals %>% filter(str_detect(id, "independent"))


# prediction of conditionals
MAP.conditionals <- predictions.MAP %>% mutate(response = utterance) %>% 
  chunk_utterances() %>% rename(utt_type = utterance, utterance = response) %>% 
  group_by(id, speaker_model, utt_type) %>% 
  filter(speaker_model == "pragmatic" & utt_type == "conditional") %>% 
  mutate(relation = case_when(str_detect(id, "if1") ~ "if1", 
                              str_detect(id, "if2") ~ "if2", 
                              str_detect(id, "independent") ~ "independent"))

MAP.conditionals %>% 
  summarize(p_hat = sum(p_hat), .groups = "drop") %>% 
  arrange(p_hat)

MAP.conditionals %>% arrange(desc(p_hat)) %>%
  group_by(relation, utterance) %>% 
  summarize(p_hat = sum(p_hat), .groups = "drop_last") %>% 
  arrange(p_hat) %>% 
  filter(relation == "if1")

freq_utts.ind_ul = df.observed %>% filter(estimate > 0.05 & id == "ind:UL")



# Posterior Predictives all speakers --------------------------------------
# prediction of conditionals
pp.speakers <- get_data_all_speakers(config_weights_relations,
                                     config_cns,
                                     extra_packages,
                                     "posterior-predictive.rds", 
                                     params$speaker_subfolder) 
hdis.pp.speakers <- mean_hdi(
  pp.speakers %>% group_by(utterance, id, speaker_model), p_hat
)
hdis.pp.speakers %>% 
  filter(str_detect(utterance, "might") & !str_detect(speaker_model, "gamma")) %>% 
  filter(p_hat < 0.05)

pp.utts <- pp.speakers %>% 
  mutate(response = utterance) %>% chunk_utterances() %>% 
  rename(utt_type = utterance, utterance = response) %>% 
  group_by(id, speaker_model, idx, utt_type) %>% 
  summarize(p_hat = sum(p_hat), .groups = "drop_last")

hdis.pragmatic <- mean_hdi(
  pp.utts %>% group_by(id, speaker_model, utt_type), p_hat
) %>% 
  filter(speaker_model %in% c("pragmatic", "pragmatic.gamma")) %>% 
  mutate(label_id = map_chr(id, get_str_contexts))
save_data(hdis.pragmatic, paste(params$mcmc_folder, 
                                "hdis-pragmatic-speakers.rds", 
                                sep = FS))

p.hdis.pragmatic <- hdis.pragmatic %>% 
  ggplot(aes(y=label_id, x = p_hat, color = speaker_model)) +
  geom_point(aes(shape = utt_type), size=4) + 
  geom_errorbarh(data = hdis.pragmatic, aes(xmin = .lower, xmax = .upper)) +
  scale_color_manual(name = "speaker model",
                     values = SPEAKER_COLORS[c("pragmatic","pragmatic.gamma")]) +
  scale_shape_manual(name = "utterance type", 
                     values = c(`might + literal` = 16, 
                                `conditional` = 17, 
                                `literal` = 15, 
                                `conjunction` = 18)) +
  scale_y_discrete(labels = label_parse()) +
  guides(shape = guide_legend(title.position = "top"), 
         color = guide_legend(title.position = "top")) +
  labs(y = "context", x = "Posterior predictive estimates with 95% HDI")
p.hdis.pragmatic  
ggsave(paste(params$mcmc_folder, "hdis-pp-pragmatic-models.png", sep=FS), 
       width = 18, plot = p.hdis.pragmatic)

# with observed estimates
# add up utterance types
df.obs <- df.observed %>% dplyr::select(estimate, id, trial, utt_type) %>% 
  group_by(utt_type, trial) %>% 
  summarize(estimate = sum(estimate), .groups = "drop_last") %>% 
  mutate(label_id = map_chr(trial, get_str_contexts)) %>% 
  add_column(speaker_model = "observed")

p.hdis.pragmatic <- hdis.pragmatic %>% 
  ggplot(aes(y=label_id, x = p_hat, color = speaker_model)) +
  geom_point(aes(shape = utt_type), size=4) + 
  geom_errorbarh(data = hdis.pragmatic, aes(xmin = .lower, xmax = .upper)) +
  geom_point(data = df.obs, aes(x=estimate, y=label_id, shape = utt_type), size = 3) + 
  scale_color_manual(name = "speaker model",
                     values = SPEAKER_COLORS[c("pragmatic","pragmatic.gamma", 
                                               "observed")]) +
  scale_shape_manual(name = "utterance type", 
                     values = c(`might + literal` = 16, 
                                `conditional` = 17, 
                                `literal` = 15, 
                                `conjunction` = 18)) +
  scale_y_discrete(labels = label_parse()) +
  guides(shape = guide_legend(title.position = "top"), 
         color = guide_legend(title.position = "top")) +
  labs(y = "context", x = "Posterior predictive estimates with 95% HDI")
p.hdis.pragmatic  
ggsave(paste(params$mcmc_folder, "hdis-pp-pragmatic-models-with-empiric.png", 
             sep=FS), width = 20, plot = p.hdis.pragmatic)


