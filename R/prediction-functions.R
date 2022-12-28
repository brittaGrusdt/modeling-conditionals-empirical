library(ModelUtils)
library(philentropy)

#' bins empirical and model tables and matches them so that for each context
#' there is is set of model tables representing participants' judged beliefs
#' Problem: not all empirical tables will be taken into account! (only those that
#' match with generated model-tables).
#' @param speaker model predictions
#' @param  tbls tibble with generated model tables, columns 'AC', 'A-C', '-AC',
#' '-A-C', 'id', 'bn_id', 'r', 'cn'
#' @param path_empiric_tbls_ids path to empirically generated probability tables
#' (csv-file)
predictions_by_matching_tables = function(speaker, tbls, path_empiric_tbls_ids){
  
  tables <- tbls %>%
    mutate(AC.round=as.integer(round(AC, 2) * 100),
           `A-C.round`=as.integer(round(`A-C`,2) * 100),
           `-AC.round`=as.integer(round(`-AC`, 2) * 100),
           `-A-C.round`=as.integer(round(`-A-C`, 2) * 100)) %>% 
    rowid_to_column("table_id")  
  
  tables.empiric_pids = match_sampled_and_empiric_tables(tables, 
                                                         path_empiric_tbls_ids) %>% 
    ungroup() %>% 
    dplyr::select(r, bn_id, empirical_id, stimulus, p_id, match.empirical)
  
  
  predictions <- left_join(speaker, tables.empiric_pids) %>%
    filter(match.empirical) %>% unnest(c(p_id)) %>% 
    dplyr::select(bn_id, utterance, probs, p_id) %>% 
    separate("p_id", into = c("subject_id", "r", "prior"), sep = "_") %>% 
    unite("id", "r", "prior", sep = "_")
  df.matches = predictions %>% dplyr::select(bn_id, subject_id, id) %>% 
    distinct() %>% group_by(subject_id, id) %>% dplyr::count() %>% arrange(desc(n))
  
  predictions.context = predictions %>% group_by(utterance, id) %>% 
    summarize(model = mean(probs), .groups = "drop_last") %>% 
    arrange(desc(model)) %>% 
    mutate(utt = utterance) %>% 
    chunk_utterances() %>% rename(utt_type = utterance, utterance = utt)

  return(list(predictions = predictions.context, matches = df.matches))
}


predictions_by_relation = function(speaker){
  predictions <- speaker %>% group_by(r, utterance) %>% 
    summarize(model = mean(probs), .groups = "drop_last") %>% 
    arrange(desc(model)) %>% 
    mutate(utt = utterance) %>% 
    chunk_utterances() %>% rename(utt_type = utterance, utterance = utt)
  return(predictions)
}

kl_wrapper = function (x, test.na = TRUE, unit = "log2", est.prob = NULL, 
          epsilon = 1e-05, mute.message = TRUE) 
{
  if (!is.matrix(x)) 
    stop("Please provide a matrix as input, e.g. with x <- rbind(vector1, vector2).", 
         call. = FALSE)
  return(distance(x = x, method = "kullback-leibler", test.na = test.na, 
                  unit = unit, est.prob = est.prob, epsilon = epsilon, 
                  mute.message = mute.message))
}



compute_kl_divergences = function(tbls.model, tbls.empiric, n_best = 10){
  model.mat <- tbls.model %>% ungroup() %>% 
    dplyr::select(AC, `A-C`, `-AC`, `-A-C`) %>% 
    as.matrix()
  kl.D = map_dfr(seq(1, nrow(tbls.empiric)), function(idx.empiric){
    print(paste("trial no", idx.empiric))
    kl.dij = map_dfr(seq(1, nrow(tbls.model)), function(idx.model){
      kl_wrapper(rbind(tbls.empiric[idx.empiric,
                                    c('AC', 'A-C', '-AC', '-A-C')] %>% as.numeric(),
                        model.mat[idx.model,])) %>% as_tibble() %>% 
        rename(kl.div = value) %>% 
        add_column(bn_id = tbls.model[idx.model, ]$bn_id)
    }) %>% arrange(kl.div) %>% 
      add_column(prolific_id = tbls.empiric[idx.empiric,]$prolific_id,
                 id = tbls.empiric[idx.empiric,]$id)
    return(kl.dij[1:n_best,])
  }) %>% 
    mutate(temp = bn_id) %>% 
    separate(temp, into = c("r.match", "prob.match", 
                            "AC.match", "A-C.match", 
                            "-AC.match", "-A-C.match"), sep = "_")
  
  n_trials = kl.D %>% group_by(prolific_id, id) %>% n_groups()
  kl.D$idx = rep(seq(1, n_best), n_trials)
  kl.D 
  return(kl.D)
}



