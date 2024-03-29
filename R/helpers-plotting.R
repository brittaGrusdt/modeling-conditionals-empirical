library(ggpubr)
library(ggthemes)
names_data <- c(`behavioral` = "darkgreen", `model` = "firebrick")

# model utterances --------------------------------------------------------
utts.model.ifs = c("A > C", "A > -C", "-A > C", "-A > -C", 
                   "C > A", "C > -A", "-C > A", "-C > -A")
utts.model.literals = c("A", "-A", "C", "-C")
utts.model.mights = c("might A", "might -A", "might C", "might -C")
utts.model.conjs = c("C and A", "C and -A", "-C and A", "-C and -A")


# mapping utterances to colors --------------------------------------------
UTT_COLORS <- c(`both blocks fall` = "darkgoldenrod1", 
                `blue falls but green does not fall` = "brown1", 
                `green falls but blue does not fall` = "brown3",
                `neither block falls` = "darkred", 
                `blue falls` = "skyblue1",
                `green falls` = "green",
                `blue does not fall` = "blue", 
                `green does not fall` = "darkgreen",
                `if blue falls green falls` = "orchid1",
                `if green falls blue falls` = "mediumvioletred",
                `if blue does not fall green does not fall` = "mediumpurple3",
                `if green does not fall blue does not fall` = "mediumorchid", 
                `if blue falls green does not fall` = "gray15",
                `if green falls blue does not fall` = "gray30",
                `if blue does not fall green falls` = "gray50",
                `if green does not fall blue falls` = "gray80", 
                `blue might fall` = "cyan",
                `green might fall` = "yellow2",
                `blue might not fall` = "turquoise4",
                `green might not fall` = "yellow4"
)


# mapping speaker models to colors ----------------------------------------
SPEAKER_COLORS <- c(
  `random` = "gray47",
  `literal` = "blue",
  `pragmatic` = "red",
  `literal.gamma` = "skyblue",
  `pragmatic.gamma` = "indianred4",
  `observed` = "turquoise4"
)

# context names -----------------------------------------------------------
get_dep_context_expression = function(trial_id){
  expr <- switch(trial_id,
                "if1_hh" = c(expression("if"[1]*":HI"), "high"),
                "if1_uh" = c(expression(paste(`if`[1], ":UI")),"unc"),
                "if1_u-Lh" = c(expression("if"[1]*":U"^-{}*"I"), "uncl"),
                "if1_lh" = c(expression("if"[1]*":LI"), "low"),
                "if2_hl" = c(expression("if"[2]*":HL"), "high"),
                "if2_ul" = c(expression("if"[2]*":UL"), "unc"),
                "if2_u-Ll" = c(expression("if"[2]*":U"^-{}*"L"), "uncl"),
                "if2_ll" = c(expression("if"[2]*":LL"), "low"))
  return(expr)
}

names.ind_contexts = list("independent_hh" = "ind:HH", 
                          "independent_ll" = "ind:LL", 
                          "independent_uh" = "ind:UH", 
                          "independent_ul" = "ind:UL", 
                          "independent_hl" = "ind:HL")

get_name_context = function(trial_id){
  if(str_detect(trial_id, "independent")){
    return(names.ind_contexts[[trial_id]])
  } else{
    return(get_dep_context_expression(trial_id))
  }
}

labels_dep_contexts <- c(
  `if1_hh` = parse(text = expression("if"[1]*":HI")),
  `if1_lh` = parse(text = expression("if"[1]*":LI")),
  `if1_u-Lh` = parse(text = expression("if"[1]*":U"^-{}*"I")),
  `if1_uh` = parse(text = expression("if"[1]*":UI")),
  `if2_hl` = parse(text = expression("if"[2]*":HL")),
  `if2_ll` = parse(text = expression("if"[2]*":LL")),
  `if2_u-Ll` = parse(text = expression("if"[2]*":U"^-{}*"L")),
  `if2_ul` = parse(text = expression("if"[2]*":UL"))
) 

# for using label_parsed!
names.contexts = list(
  "if1_uh" = "'if'[1]*':UI'",
  "if1_hh" = "'if'[1]*':HI'",
  "if1_lh" = "'if'[1]*':LI'",
  "if1_u-Lh" = "'if'[1]*':U'^-{}*'I'",
  "if2_u-Ll" = "'if'[2]*':U'^-{}*'L'",
  "if2_ul" = "'if'[2]*':UL'",
  "if2_hl" = "'if'[2]*':HL'",
  "if2_ll" = "'if'[2]*':LL'",
  "independent_hh" = "ind:HH", 
  "independent_ll" = "ind:LL", 
  "independent_uh" = "ind:UH", 
  "independent_ul" = "ind:UL", 
  "independent_hl" = "ind:HL"
)
get_str_contexts = function(trial){
  return(names.contexts[[trial]])
}

# functions ---------------------------------------------------------------
plot_correlation = function(results.joint, ncol, label.x = NA, label.y = 0.6){
  p_scatter = results.joint %>% 
    ggscatter(y = "behavioral", x = "model", 
              add = "reg.line", conf.int = TRUE, cor.method = "pearson", 
              ylab = "Data", xlab = "Predictions") +
    geom_point(size=2, 
               aes_string(y="behavioral", x="model", color = "utterance", 
                          fill = "utterance")
               ) +
    guides(fill = guide_legend(title.position = "top"), 
           color = guide_legend(title.position = "top")) +
    facet_wrap(~label_id, labeller = label_parsed, ncol = ncol) +
    scale_color_manual(name = "utterance", values = UTT_COLORS) +
    scale_fill_manual(name = "utterance", values = UTT_COLORS) +
    theme_clean(base_size = 26) + 
    theme(legend.position = "top", text = element_text(size = 26))
  if(!is.na(label.x)) {
    p_scatter <- p_scatter + 
      stat_cor(label.x = label.x, label.y = label.y, size = 6)
  } else {
    p_scatter <- p_scatter + stat_cor(label.y = label.y, size = 6)
  }
  return(p_scatter)
}


plot_model_vs_data_bars = function(df.joint, tit = "", by_utt_type=F){
  df.long <- df.joint %>%
    pivot_longer(cols = c("model", "behavioral"), names_to = "data",
                 values_to = "probability") %>% 
    arrange(data, id, desc(probability)) %>% 
    group_by(data, id) %>% mutate(cdf = cumsum(probability))
  if(by_utt_type){
    df.long <- df.long %>% mutate(utt=utterance) %>% chunk_utterances() %>% 
      rename(utt_type = utterance, utterance = utt) %>% 
      group_by(id, relation, data, utt_type) %>% 
      summarize(probability = sum(probability), .groups = "drop_last") %>% 
      rename(utterance = utt_type)
  }
  plots = map(df.long$relation %>% unique(), function(rel){
    p <- df.long %>% filter(relation == !!rel) %>% 
      filter(cdf <= 0.99) %>% 
      ggplot(aes(y=utterance, x=probability, fill=data)) + 
      geom_bar(stat="identity", position = position_dodge()) +
      labs(x = "condition", title = tit) + 
      facet_wrap(~id, scales = "free_y") +
      scale_fill_manual(values = names_data)
    return(p)
  })
  return(plots)
}

# returns matrix nb.drawn samples x N=nb.participants, 
# for a single sampled value, given in arg 'col'
format_sampled_ps <- function(df, col) {
  df %>% dplyr::select(idx_sample, id_subj, !!col) %>% 
    group_by(idx_sample) %>% 
    pivot_wider(names_from="id_subj", values_from=col) %>% 
    ungroup() %>% dplyr::select(-idx_sample) %>% as.matrix()
}

