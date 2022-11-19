plot_correlation = function(results.joint, color_utt = F){
  if(color_utt) {
    p_scatter = results.joint %>% 
      ggscatter(x = "behavioral", y = "model", color = "utterance", 
                add = "reg.line", conf.int = TRUE, cor.coef = TRUE, 
                cor.method = "pearson", xlab = "Empirical observations",
                ylab = "Model predictions") +
        geom_point(size=1.5, aes_string(x="behavioral", y="model", 
                                        color="utterance")) + 
        theme(legend.position = "none")
  } else {
    p_scatter = results.joint %>% 
      ggscatter(x = "behavioral", y = "model", add = "reg.line",
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
                xlab = "Empirical observations", ylab = "Model predictions") +
      geom_point(size=1.5, aes_string(x="behavioral", y="model", color = "utterance"))
  }

  return(p_scatter)
}