library(here)
library(boot)
library(ggforce)
library(xtable)
library(latex2exp)
library(scales)
library(tidyverse)
library(ExpDataWrangling)
library(ggthemes)
theme_set(theme_clean(base_size=20) + theme(legend.position = "top"))

path_cleaned_data = here("data", "cleaned-data.csv")

params <- config::get()
plot_dir = here("data", "figs")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = T)


# Data --------------------------------------------------------------------
cleaned_data = read_csv(path_cleaned_data) %>% get_controlled_factors()

# select smoothed or original slider ratings!
pe_data.long =  cleaned_data %>%
  dplyr::select(prolific_id, id, relation, pe_task, slider,
                prior_blue, prior_green) %>% 
  filter(!is.na(slider)) %>% 
  group_by(id, prolific_id) %>% 
  mutate(slider = case_when(slider == "bg" ~ "ac", 
                            slider == "b" ~ "a-c",
                            slider == "g" ~ "-ac",
                            slider == "none" ~ "-a-c")) %>% 
  rename(world = slider) %>% 
  mutate(world = factor(world, levels = c("ac", "a-c", "-ac", "-a-c")), 
         id = str_replace(id, "independent", "ind"),
         id = factor(id, levels = c("if1_hh", "if1_uh", "if1_u-Lh",
                                    "if1_lh", "if2_hl", "if2_ul", "if2_u-Ll",
                                    "if2_ll", "ind_hh", "ind_uh", "ind_hl",
                                    "ind_ul", "ind_ll")))

pe_data.wide = pe_data.long %>% 
  pivot_wider(names_from = "world", values_from = "pe_task") %>% 
  mutate(a = ac + `a-c`, c = ac + `-ac`)

data.uc = cleaned_data %>%  
  mutate(utterance = utt.standardized) %>% 
  chunk_utterances() %>%
  rename(utt_type = utterance) %>%  #, utterance = utt.standardized) %>% 
  mutate(utt_type = as.character(utt_type)) %>% 
  filter(human_exp2 == 1) %>% 
  dplyr::select(prolific_id, id, relation, relation_type, prior_blue, prior_green,
                utt_type, uc_task, utt.standardized, pe_task.smooth, pe_task, 
                cost.uc) %>% 
  mutate(conditional = case_when(utt_type != "conditional" ~ 0, T ~ 1),
         prior_blue_agg = case_when(as.character(prior_blue) == "unc" ~ "unc", 
                                    as.character(prior_blue) == "uncl" ~ "uncl",
                                    T ~ "confident"),
         prior_blue_agg = factor(prior_blue_agg, levels = c("unc", "uncl", 
                                                            "confident"),
                                 labels = c("unc", "uncl", "confident")),
         prior_blue = factor(prior_blue,
                             levels = c("low", "uncl", "unc", "high"),
                             labels = c("L", "U^-", "U", "H"))
  ) %>% 
  rename(subject_id = prolific_id)

# Bootstrap PE-task data ---------------------------------------------------
get_mean_estimates = function(data, i){
  d2 <- data[i,] %>% as.matrix()
  return(colMeans(d2))
}
N = 1000

mat <- pe_data.wide %>% group_by(id, relation) %>% 
  rename(bg=ac, b=`a-c`, g=`-ac`, none=`-a-c`) %>% 
  ungroup() %>% dplyr::select(bg, b, g, none) %>% as.matrix()
pe_data.wide$y <- mat

bootstrapped.pe = group_map(pe_data.wide %>% 
                                group_by(id, relation, prior_blue, prior_green), 
                              function(df.grp, grp){
    set.seed(1234)
    samples <- boot(df.grp$y %>% as_tibble(),
                    statistic = get_mean_estimates, R=N)
    bg = boot.ci(samples, index=c(1), type = "basic", conf=c(0.95))$basic[4:5]
    b = boot.ci(samples, index=c(2), type = "basic", conf=c(0.95))$basic[4:5]
    g = boot.ci(samples, index=c(3), type = "basic", conf=c(0.95))$basic[4:5]
    none = boot.ci(samples, index=c(4), type = "basic", conf=c(0.95))$basic[4:5]
    
    cis = rbind(bg, b, g, none)
    colnames(cis) <- c("lower", "upper")
    cis %>% as_tibble(rownames=NA) %>% rownames_to_column("world") %>% 
      add_column(pblue = grp$prior_blue, relation = grp$relation, 
                 pgreen = grp$prior_green, id=grp$id)
  }) %>% bind_rows() %>%  add_column(data = "empiric") %>% 
  mutate(world = factor(world, levels = c("bg", "b", "g", "none"),
                        labels = c("bg", "b¬g", "¬bg", "¬b¬g"))) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent"), 
                           labels = c(parse(text = expression("if"[1])), 
                                      parse(text = expression("if"[2])), 
                                      "ind")))

pe_task.means =  pe_data.long %>% group_by(id, relation, world) %>% 
  summarize(pe_task = mean(pe_task), .groups = "drop_last") %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))
pe_task.means$world <- recode_factor(pe_task.means$world, 
                                     `ac` = "bg", `a-c` = "b¬g", 
                                     `-ac` = "¬bg", `-a-c` = "¬b¬g")

p.means_pe = pe_task.means %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent"), 
                           labels = c(parse(text = expression("if"[1])), 
                                      parse(text = expression("if"[2])), 
                                      "ind")
                           )) %>% 
  ggplot(aes(x = id,  color = world, group = world)) +
  geom_point(aes(y = pe_task)) + geom_line(aes(y = pe_task)) +
  scale_x_discrete(labels = c(`if1_hh` = parse(text = TeX('$HI$')),
                       `if1_lh` = parse(text = TeX('$LI$')),
                       `if1_u-Lh` = parse(text = TeX('$U^-I$')),
                       `if1_uh` = parse(text = TeX('$UI$')),
                       `if2_hl` = parse(text = TeX('$HL$')),
                       `if2_ll` = parse(text = TeX('$LL$')),
                       `if2_u-Ll` = parse(text = TeX('$U^-L$')),
                       `if2_ul` = parse(text = TeX('$UL$')),
                       `ind_hh` = 'HH',
                       `ind_uh` = 'UH',
                       `ind_hl` = 'HL',
                       `ind_ul` = 'UL',
                       `ind_ll` = 'LL'
                     )
  ) +
  geom_errorbar(data = bootstrapped.pe, aes(ymin = lower, ymax = upper)) +
  facet_wrap(relation ~ world, scales = "free_x", ncol = 4,  
             labeller = labeller(relation = label_parsed)) +
  labs(y = "Mean estimates PE-task", x = "condition (context)")
p.means_pe
ggsave(paste(plot_dir, "mean-estimates-pe-task.png", sep = FS), 
       p.means_pe, width = 14, height = 11)

p.means_pe.ribbons = pe_task.means %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent"), 
                           labels = c(parse(text = expression("if"[1])), 
                                      parse(text = expression("if"[2])), 
                                      "ind"))) %>% 
  ggplot(aes(x = id,  color = world, group = world)) +
  geom_point(aes(y = pe_task)) + geom_line(aes(y = pe_task)) +
  scale_x_discrete(labels = 
                     c(`if1_hh` = parse(text = TeX('$HI$')),
                       `if1_lh` = parse(text = TeX('$LI$')),
                       `if1_u-Lh` = parse(text = TeX('$U^-I$')),
                       `if1_uh` = parse(text = TeX('$UI$')),
                       `if2_hl` = parse(text = TeX('$HL$')),
                       `if2_ll` = parse(text = TeX('$LL$')),
                       `if2_u-Ll` = parse(text = TeX('$U^-L$')),
                       `if2_ul` = parse(text = TeX('$UL$')),
                       `ind_hh` = 'HH',
                       `ind_uh` = 'UH',
                       `ind_hl` = 'HL',
                       `ind_ul` = 'UL',
                       `ind_ll` = 'LL'
                     )) +
  geom_errorbar(data = bootstrapped.pe, aes(ymin = lower, ymax = upper)) +
  geom_ribbon(data = bootstrapped.pe, 
              aes(ymin = lower, ymax = upper, fill=world),
              alpha=0.5) +
  facet_wrap(relation ~ world, scales = "free_x", ncol = 4, 
             labeller = labeller(relation = label_parsed)) +
  labs(y = "Mean estimates PE-task", x = "condition (context)") +
  theme(text = element_text(size=18))
p.means_pe.ribbons
ggsave(paste(plot_dir, "mean-estimates-pe-task-ribbons.png", sep = FS), 
       p.means_pe.ribbons, width = 16, height = 11)



# UC-task data ------------------------------------------------------------

# 1. standardized utterances selection frequencies
df.utt_freq = data.uc %>% 
  dplyr::select(subject_id, id, utt.standardized, utt_type) %>% 
  group_by(utt.standardized, utt_type) %>% dplyr::count() %>% arrange(n)
df.utt_freq$utt_f = factor(df.utt_freq$utt.standardized, 
                           levels = df.utt_freq$utt.standardized)

p.utt_freq = df.utt_freq %>% ggplot(aes(y = utt_f, x = n)) + 
  geom_bar(aes(fill = utt_type), stat = "identity") +
  scale_fill_brewer(name = "utterance type", palette = "Set2") +
  geom_text(aes(label = n), hjust = -0.2, size=6) +
  labs(y = "standardized utterance", x = "number selections") #+
# theme(text = element_text(size = 18))
p.utt_freq
ggsave(paste(plot_dir, "utterance-frequencies.png", sep = FS), p.utt_freq, 
       width = 14, height = 6)

# 2. selections of standardized conditionals
df.if = data.uc %>% filter(utt_type == "conditional") %>% 
  group_by(id, utt.standardized, uc_task) %>% 
  dplyr::count() %>% arrange(desc(n))

df.if.standardized =  df.if %>% 
  group_by(id, utt.standardized) %>% 
  summarize(n = sum(n), .groups = "drop") %>% arrange(desc(n))
positions_ids = df.if.standardized$id %>% unique()
p.ifs.standardized = df.if.standardized %>% 
  ggplot(aes(y = n, x = id, fill = utt.standardized)) +
  geom_bar(stat = "identity") +
  labs(x = "condition (context)", y = "count") + 
  scale_x_discrete(limits = positions_ids, 
                   labels = 
                     c(`if1_hh` = parse(text = expression("if"[1]*":HI")),
                       `if1_lh` = parse(text = expression("if"[1]*":LI")),
                       `if1_u-Lh` = parse(text = expression("if"[1]*":U"^-{}*"I")),
                       `if1_uh` = parse(text = expression("if"[1]*":UI")),
                       `if2_hl` = parse(text = expression("if"[2]*":HL")),
                       `if2_ll` = parse(text = expression("if"[2]*":LL")),
                       `if2_u-Ll` = parse(text = expression("if"[2]*":U"^-{}*"L")),
                       `if2_ul` = parse(text = expression("if"[2]*":UL")),
                       `independent_hh` = 'ind:HH',
                       `independent_uh` = 'ind:UH',
                       `independent_hl` = 'ind:HL',
                       `independent_ul` = 'ind:UL',
                       `independent_ll` = 'ind:LL'
                     )) +
  scale_fill_brewer(name = "selected (standardized) conditional", 
                    palette = "Set2") +
  guides(fill = guide_legend(ncol=2, title.position = "top"))# +
p.ifs.standardized
ggsave(paste(plot_dir, "conditionals-uc-task.png", sep=FS), p.ifs.standardized,
       width = 11, height = 6)



# 3. Ratio of conditionals dependent on relation and prior conditions
get_mean_estimates = function(data, i){
  d2 <- data[i,]
  return(mean(d2$conditional))
}

N = 1000
pconditionals_bootstrapped = group_map(
  data.uc %>% group_by(relation, prior_blue), function(df.grp, grp){
    samples <- boot(df.grp, 
                    statistic = get_mean_estimates, R=N)$t %>% sort()
    samples <- boot(df.grp, statistic = get_mean_estimates, R=N)
    cis <- boot.ci(samples, type = "basic", conf=c(0.95))$basic[4:5]
    
    tibble(lower=cis[1], upper=cis[2]) %>% 
      add_column(prior_blue = grp$prior_blue, relation=grp$relation)
  }) %>% bind_rows()  


# Model relation:dependent vs. independent, confident vs. uncertain
# plot data 
pd <- position_dodge(0.2)
p.ratio_conditionals <- data.uc %>%
  group_by(relation, prior_blue) %>% 
  summarize(N = n(), num_conditionals = sum(conditional), .groups = "drop_last") %>% 
  mutate(ratio = num_conditionals / N) %>% 
  ggplot(aes(x = relation, group = prior_blue, color = prior_blue)) +
  geom_point(position = pd, aes(y=ratio)) + 
  geom_line(position = pd, aes(y=ratio)) +
  geom_errorbar(data = pconditionals_bootstrapped,
                position = pd, aes(ymin = lower, ymax = upper)
  ) +
  scale_x_discrete(limits = c("if1", "if2", "independent"), 
                   labels = 
                     c(`if1` = parse(text = expression("if"[1])),
                       `if2` = parse(text = expression("if"[2])),
                       `independent` = 'independent'
                     )) +
  scale_color_discrete(labels = unname(TeX(c("L", "$U^-$", "U", "H")))) +
  ylab("ratio selected conditionals")
p.ratio_conditionals
ggsave(paste(plot_dir, "ratio_conditionals_empirical.png", sep=FS), 
       plot = p.ratio_conditionals)


# Joint data PE- and UC-task ----------------------------------------------
pe_uc_data = left_join(
  pe_data.wide %>% dplyr::select(-y),
  data.uc %>% rename(prolific_id = subject_id) %>% 
    dplyr::select(prolific_id, id, relation, utt_type, uc_task, utt.standardized) %>% 
    mutate(id = str_replace(id, "independent", "ind"))
) %>% rename(bg = ac, b = `a-c`, g = `-ac`, none = `-a-c`) %>% 
  add_probs() %>% rename(ac = bg, `a-c` = `b`, `-ac` = g, `-a-c` = none) %>% 
  compute_utt_probs() %>% dplyr::select(-starts_with("p_")) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))

df.means_pe_utts.by_rel = pe_uc_data %>% 
  group_by(utt.standardized, utt_type, relation) %>%
  summarize(mean_pe = mean(pe_selected_utt), .groups = "drop_last") %>% 
  arrange(desc(mean_pe)) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))
# conjunctions
df.means_pe_utts.by_rel %>% filter(utt_type == "conjunction") %>% arrange(mean_pe)
# overall mean / median by relation
df.means_pe_utts.by_rel %>% filter(utt_type == "conjunction") %>% 
  group_by(utt.standardized) %>% 
  summarize(mean=mean(mean_pe), median = median(mean_pe))

# selections of utterances
df.means_pe_utts.by_rel %>% filter(utt.standardized == "if blue falls green falls")
df.means_pe_utts.by_rel %>% filter(utt.standardized == "blue might fall" & relation == "if2")
df.means_pe_utts.by_rel %>% filter(utt.standardized == standardized.sentences$bg)
df.means_pe_utts.by_rel %>% filter(utt.standardized == standardized.sentences$none)
df.means_pe_utts.by_rel %>% filter(utt.standardized == standardized.sentences$b)
df.means_pe_utts.by_rel %>% filter(utt.standardized == standardized.sentences$g)

# min average value for literals except 'green falls' in if2 selected only once
df.means_pe_utts.by_rel %>% filter(utt_type == "literal") %>% 
  filter(!(utt.standardized == "green falls" & relation == "if2")) %>% 
  pull(mean_pe) %>% min()

# mean estimate in pe-task for a single selected utterance
get_mean_estimates = function(data, i){
  d2 <- data[i,]
  return(mean(d2$pe_selected_utt))
}
N = 1000
data_bootstrapped_pe_uc = group_map(
  pe_uc_data %>% group_by(relation, utt.standardized), function(df.grp, grp){
    samples <- boot(df.grp, statistic = get_mean_estimates, R=N)$t %>% sort()
    tibble(x_hat = list(samples), 
           lower = samples[round(0.025*N)], 
           upper = samples[round(0.975*N)],
           relation = grp$relation, utt.standardized = grp$utt.standardized)
  }) %>% bind_rows() %>% mutate(utterance = utt.standardized) %>% 
  chunk_utterances() %>% rename(utt_type = utterance) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))

df.utt_freq_by_rel = data.uc %>% ungroup() %>% 
  dplyr::select(subject_id, relation, utt.standardized, utt_type) %>% 
  group_by(relation, utt.standardized, utt_type) %>% dplyr::count() %>% 
  arrange(n) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))

p.pe_means_pe_utts = df.means_pe_utts.by_rel %>% 
  ggplot(aes(y = utt.standardized, color = relation)) + 
  geom_point(aes(x = mean_pe),
             position = position_dodge(width=0.5)) + 
  facet_wrap(~utt_type, scales = "free_y") +
  geom_errorbar(data = data_bootstrapped_pe_uc, 
                aes(xmin = lower, xmax = upper), 
                position = position_dodge(width = 1)) +
  geom_text(data = df.utt_freq_by_rel, aes(x = 0.1, label = n), 
            position = position_dodge(width=1), size=6) +
  scale_color_brewer(name = "relation", palette = "Set2", 
                     labels = c(`if1` = parse(text = expression("if"[1])),
                                `if2` = parse(text = expression("if"[2])),
                                `independent` = 'independent')) +
  labs(x = "Mean estimate of probability for described event", 
       y = "selected (standardized) utterance")
p.pe_means_pe_utts

ggsave(paste(plot_dir, "pe_means_pe_utts.png", sep = FS), p.pe_means_pe_utts,
       width = 14, height = 8)

# bootstrapped utterance ratios for each context
get_mean_uc_ratios = function(data, i, utt){
  d2 <- data[i,]
  return(mean(d2$utt.standardized == utt))
}
N = 1000
bootstrapped_uc_ratios_ci = group_map(
  pe_uc_data %>% mutate(id = str_replace(id, "ind", "independent")) %>% 
    dplyr::select(prolific_id, id, utt.standardized) %>%  group_by(id), 
  function(df.id, grp.id){
    message(grp.id$id)
    df <- group_map(df.id %>% group_by(utt.standardized), function(df.utt, grp.utt){
      samples <- boot(df.id, statistic = get_mean_uc_ratios, R=N, 
                      utt=grp.utt$utt.standardized)$t %>% sort()
      tibble(x_hat = list(samples), 
             estimate = mean(samples),
             lower = samples[round(0.025*N)], 
             upper = samples[round(0.975*N)],
             id = grp.id$id, utt.standardized = grp.utt$utt.standardized)
    }) %>% bind_rows()
    return(df)
  }) %>% bind_rows() 
save_data(bootstrapped_uc_ratios_ci, here(params$dir_data, 
                                        "bootstrapped_uc_ratios_ci.rds"))

# estimated P(blue) vs. P(green) when conditional selected ----------------
df.utt_freq_by_rel = data.uc %>% ungroup() %>% 
  dplyr::select(subject_id, relation, utt.standardized, utt_type) %>% 
  group_by(relation, utt.standardized, utt_type) %>% dplyr::count() %>% 
  arrange(n) %>% 
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent")))

to_remove = df.utt_freq_by_rel %>% 
  filter(utt.standardized == standardized.sentences$if_bg & n <= 2)

df.ifbg = pe_uc_data %>% 
  mutate(na = `-ac` + `-a-c`, nc = `a-c` + `-a-c`) %>% 
  dplyr::select(prolific_id, id, relation, a, c, na, nc,
                uc_task, utt_type, utt.standardized) %>%
  mutate(
    p_ant = case_when(startsWith(utt.standardized, "if blue falls") |
                      startsWith(utt.standardized, "if blue does not fall") ~ a,
                      startsWith(utt.standardized, "if green falls") |
                      startsWith(utt.standardized, "if green does not fall") ~ c),
    p_cons = case_when(endsWith(utt.standardized, "blue falls") |
                       endsWith(utt.standardized, "blue does not fall") ~ a,
                       endsWith(utt.standardized, "green falls") |
                       endsWith(utt.standardized, "green does not fall") ~ c)) %>%
  filter(utt.standardized == standardized.sentences$if_bg) 

p.ant_vs_cons_ifbg = anti_join(df.ifbg, to_remove) %>%  
  mutate(relation = factor(relation, levels = c("if1", "if2", "independent"), 
                           labels = c(parse(text = expression("if"[1])), 
                                      parse(text = expression("if"[2])), 
                                      "ind"))) %>% 
  ggplot(aes(x = p_ant, y = p_cons)) + 
  geom_density_2d_filled(contour_var = "ndensity") +
  geom_point(alpha = 0.5) + xlim(0, 1) + ylim(0, 1) +
  facet_wrap(relation~utt.standardized, labeller = labeller(relation = label_parsed)) +
  labs(x = "Estimate P(blue falls)", y = "Estimate P(green falls)")
p.ant_vs_cons_ifbg

ggsave(paste(plot_dir, "pant_vs_pcons.png", sep = FS), p.ant_vs_cons_ifbg, 
       width = 11, height = 7)




# Utterance cost ----------------------------------------------------------
# utterances: both fall & neither falls
df.cost_bg_nbng = data.uc %>%
  mutate(cost.uc = factor(cost.uc, levels = seq(2, 7))) %>% 
  filter(utt.standardized %in% 
           c(standardized.sentences$bg, standardized.sentences$none))

df.cost_bg_nbng.sum = df.cost_bg_nbng %>% 
  group_by(utt.standardized, cost.uc, uc_task) %>% 
  summarize(n = n(), .groups = "drop_last") %>% arrange(desc(n))
df.cost_bg_nbng.sum$utt_f = factor(df.cost_bg_nbng.sum$utt.standardized, 
                                   levels = c(df.cost_bg_nbng.sum$utt.standardized %>% unique()))

p.bg_nbng_cost = df.cost_bg_nbng.sum %>% 
  ggplot(aes(x=n, y = utt_f, fill = cost.uc)) + 
  geom_bar(stat = "identity") + 
  scale_fill_brewer(name = "utterance cost (nb of clicks)", palette = 'OrRd') +
  labs(x = "number selections", y = "standardized utterance") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(text = element_text(size=18))
p.bg_nbng_cost


ggsave(paste(plot_dir, "bg_nbng_cost.png", sep = FS), p.bg_nbng_cost, 
       width = 12, height = 3)






