get_assertable_utterances = function(tbls.wide){
  tbls.wide %>% 
    rename(bg = AC, b = `A-C`, g = `-AC`, none = `-A-C`) %>% 
    add_probs() %>% 
    rename(p_AC = bg, `p_A-C` = b, `p_-AC` = g, `p_-A-C` = none) %>% 
    pivot_longer(starts_with("p_"), names_to = "utt", values_to = "prob") %>% 
    mutate(utt = case_when(utt == "p_AC" ~ standardized.sentences$bg, 
                           utt == "p_A-C" ~ standardized.sentences$b, 
                           utt == "p_-AC" ~ standardized.sentences$g,
                           utt == "p_-A-C" ~ standardized.sentences$none, 
                           utt == "p_a" ~ standardized.sentences$only_b,
                           utt == "p_c" ~ standardized.sentences$only_g,
                           utt == "p_na" ~ standardized.sentences$only_nb,
                           utt == "p_nc" ~ standardized.sentences$only_ng, 
                           utt == "p_c_given_a" ~ standardized.sentences$if_bg,
                           utt == "p_c_given_na" ~ standardized.sentences$if_nbg,
                           utt == "p_nc_given_a" ~ standardized.sentences$if_bng,
                           utt == "p_nc_given_na" ~ standardized.sentences$if_nbng,
                           utt == "p_a_given_c" ~ standardized.sentences$if_gb,
                           utt == "p_a_given_nc" ~ standardized.sentences$if_ngb,
                           utt == "p_na_given_c" ~ standardized.sentences$if_gnb,
                           utt == "p_na_given_nc" ~ standardized.sentences$if_ngnb, 
                           utt == "p_likely_a" ~ standardized.sentences$might_b,
                           utt == "p_likely_c" ~ standardized.sentences$might_g,
                           utt == "p_likely_na" ~ standardized.sentences$might_nb,
                           utt == "p_likely_nc" ~ standardized.sentences$might_ng
    ))
}