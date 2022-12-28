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
      stat_cor(aes(color = utterance), label.x = label.x, label.y = label.y)
  } else {
    p_scatter <- p_scatter + 
      stat_cor(aes(color = utterance), label.x = label.x)
  }

  return(p_scatter)
}