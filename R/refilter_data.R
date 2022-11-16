library(here)
library(tidyverse)
library(ExpDataWrangling)
library(ModelUtils)


# Data --------------------------------------------------------------------
params <- config::get()
path_cleaned_data = here(params$dir_data, params$fn_cleaned_data)
cleaned_data = read_csv(path_cleaned_data)

pe_data.long =  cleaned_data %>%
  dplyr::select(prolific_id, id, pe_task, slider) %>% 
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

uc_data = cleaned_data %>% filter(!is.na(uc_task)) %>% 
  dplyr::select(prolific_id, id, uc_task, utt.standardized, 
                custom_response, cost.uc, RT.uc_task, starts_with("intended")) %>%
  group_by(prolific_id, id) %>% 
  mutate(utterance = uc_task, 
         id = str_replace(id, "independent", "ind"), 
         id = factor(id, levels = c("if1_hh", "if1_uh", "if1_u-Lh",
                                    "if1_lh", "if2_hl", "if2_ul", "if2_u-Ll",
                                    "if2_ll", "ind_hh", "ind_uh", "ind_hl",
                                    "ind_ul", "ind_ll"))) %>% 
  chunk_utterances() %>% rename(utt_type = utterance)

# Participants Ratings vs. selected utterance -----------------------------
pe_uc_data = left_join(
  pe_data.wide,
  uc_data %>% dplyr::select(-custom_response, -cost.uc, -RT.uc_task)
) %>% rename(bg = ac, b = `a-c`, g = `-ac`, none = `-a-c`) %>% 
  add_probs() %>% rename(ac = bg, `a-c` = `b`, `-ac` = g, `-a-c` = none) %>% 
  compute_utt_probs() %>% dplyr::select(-starts_with("p_")) %>% 
  translate_standardized2model()

# exclude participants when an event other than the selected conjunction is rated
# as more likely
conjunctions = pe_uc_data %>%
  dplyr::select(prolific_id, id, ac, `a-c`, `-ac`, `-a-c`, 
                utt.standardized, utt_type, pe_selected_utt) %>% 
  filter(utt_type == "conjunction") %>%
  pivot_longer(cols = c(ac, `a-c`, `-ac`, `-a-c`), names_to = "world", values_to = "rating") %>% 
  mutate(p_max = max(rating)) %>% 
  pivot_wider(names_from = "world", values_from = "rating") %>%
  filter(pe_selected_utt != p_max)

# selected conditional is A->C, but A-C is rated as more likely than AC
conditionals <- pe_uc_data %>%
  filter(utt_type == "conditional") %>% 
  mutate(remove = case_when(utterance == "A > C"  ~ ac < `a-c`, 
                            utterance == "-A > -C" ~ `-a-c`< `-ac`,
                            utterance == "C > A" ~ ac < `-ac`,
                            utterance == "-C > A" ~ `a-c` < `-a-c`, 
                            utterance == "A > -C" ~ `a-c` < `ac`,
                            utterance == "-A > C" ~ `-ac` < `-a-c`
                            )) %>% 
  filter(remove)

# literals: 
# exclude participants when an event other than that of the selected literal utterance
# is rated as more likely (eg. 'blue falls' said but 'green falls' rated as more likely)
literals = pe_uc_data %>%
  dplyr::select(prolific_id, id, a, c, `a-c`, `-a-c`, `ac`, `-ac`, utt.standardized, utt_type, 
                utterance, pe_selected_utt) %>% 
  filter(utt_type == "literal") %>%
  mutate(na = `-ac` + `-a-c`, nc = `a-c` + `-a-c`) %>% 
  dplyr::select(-`ac`, -`a-c`, -`-ac`, -`-a-c`) %>% 
  pivot_longer(cols = c(a, c, na, nc), names_to = "marginal", values_to = "rating") %>% 
  mutate(p_max = max(rating)) %>% 
  pivot_wider(names_from = "marginal", values_from = "rating") %>%
  filter(pe_selected_utt != p_max)


trials_out = bind_rows(literals %>% dplyr::select(prolific_id, id), 
                       conjunctions  %>% dplyr::select(prolific_id, id),
                       conditionals %>% dplyr::select(prolific_id, id))


recleaned_data = anti_join(cleaned_data, trials_out)


# Save data ---------------------------------------------------------------
Sys.setenv(R_CONFIG_ACTIVE = "recleaned_data")
par.recleaned = config::get()
write_csv(recleaned_data, paste(params$dir_data, par.recleaned$fn_recleaned_data, 
                                sep = FS))

# for use with model save tables
tables = recleaned_data %>% filter(!is.na(slider)) %>% 
  dplyr::select(prolific_id, id, prob, pe_task.smooth) %>% 
  pivot_wider(names_from = "prob", values_from = "pe_task.smooth")
path_empiric_tbl_ids = paste(par.recleaned$dir_data,
                             par.recleaned$fn_tbls_empiric_pids, sep = FS)
df = generate_and_save_empiric_tbl_ids(tables, path_empiric_tbl_ids)






