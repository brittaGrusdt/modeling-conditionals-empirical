library(ggpubr)
names_data <- c(`behavioral` = "darkgreen", `model` = "firebrick")

plot_correlation = function(results.joint, 
                            color = "utterance", shape = "relation",
                            label.x = 0.1, label.y = NA){
  p_scatter = results.joint %>% 
    ggscatter(x = "behavioral", y = "model", color = color, shape = shape,
              add = "reg.line", conf.int = TRUE, cor.method = "pearson", 
              xlab = "Empirical observations", ylab = "Model predictions") +
      geom_point(size=1.5, aes_string(x="behavioral", y="model", 
                                      color=color, shape=shape)) +
    guides(fill = "none")
  if(!is.na(label.y)) {
    p_scatter <- p_scatter + 
      stat_cor(aes(color = color), label.x = label.x, label.y = label.y)
  } else {
    p_scatter <- p_scatter + 
      stat_cor(aes(color = color), label.x = label.x)
  }

  return(p_scatter)
}


plot_model_vs_data_bars = function(df.joint, tit = "", by_utt_type=F){
  df.long <- df.joint %>%
    pivot_longer(cols = c("model", "behavioral"), names_to = "data",
                 values_to = "probability")
  if(by_utt_type){
    df.long <- df.long %>% mutate(utt=utterance) %>% chunk_utterances() %>% 
      rename(utt_type = utterance, utterance = utt) %>% 
      group_by(id, relation, data, utt_type) %>% 
      summarize(probability = sum(probability), .groups = "drop_last") %>% 
      rename(utterance = utt_type)
  }
  plots = map(df.long$relation %>% unique(), function(rel){
    p <- df.long %>% filter(relation == !!rel) %>% 
      ggplot(aes(y=utterance, x=probability, fill=data)) + 
      geom_bar(stat="identity", position = position_dodge()) +
      labs(x = "condition", title = tit) + 
      facet_wrap(~id, scales = "free") +
      scale_fill_manual(values = names_data)
    return(p)
  })
  return(plots)
}
