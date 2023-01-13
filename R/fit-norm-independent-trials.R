library(here)
library(rwebppl)
library(tidyverse)
library(truncnorm)
library(greta)
library(ExpDataWrangling)
library(ModelUtils)
library(latex2exp)

active_config = "context_free_prior"
Sys.setenv(R_CONFIG_ACTIVE = active_config)

params <- config::get()
# Data --------------------------------------------------------------------
data.behav <- read_csv(here(params$dir_data, "cleaned-data.csv")) %>% 
  dplyr::select(prolific_id, id, utt.standardized, uc_task, pe_task.smooth, slider) %>% 
  translate_standardized2model() 
data.uc = data.behav %>%  filter(!is.na(uc_task)) %>% 
  dplyr::select(prolific_id, id, utterance) %>% 
  add_column(n = 1) %>% group_by(id, utterance)
data.pe = data.behav %>% 
  dplyr::select(prolific_id, id, utt.standardized, pe_task.smooth) %>% 
  pivot_wider(names_from = "utt.standardized", values_from = "pe_task.smooth") %>% 
  rename(AC = `both blocks fall`, 
         `A-C` = `blue falls but green does not fall`, 
         `-AC` = `green falls but blue does not fall`, 
         `-A-C` = `neither block falls`)

pe.diff_ind = data.pe %>% 
  group_by(prolific_id, id) %>% 
  transmute(pa = AC + `A-C`, pc = `AC` + `-AC`, p_ac.ind = pa * pc, 
            diff = AC - p_ac.ind) %>% 
  separate("id", into = c("relation", "prior"), sep = "_") %>% 
  mutate(abs_diff = abs(diff)) %>% ungroup() %>% 
  arrange(abs_diff) 

theme_set(theme_bw(base_size = 12))
p = pe.diff_ind %>% ggplot() + 
  geom_density(aes(x = diff, color = relation)) +
  labs(x = TeX(r"( $P(b,g) - P(b)\cdot P(g)$ )"),
       title = "Deviation from perfect stochastic independence") +
  theme(legend.position = c(0.8, 0.8), legend.direction = "horizontal")
p

# fit sd normal distribution ----------------------------------------------
draws.ind = pe.diff_ind %>% filter(relation == "independent") %>% rowid_to_column()
x = 2/3
# x% of smallest absolute values
vals.to_fit = draws.ind %>% ungroup() %>% filter(rowid <= x * n()) %>% 
  pull(diff)

sd <- greta::uniform(0, 0.05, 1) 
greta::distribution(vals.to_fit) <- greta::normal(0, sd)
model = greta::model(c(sd))

fit = opt(model)
fitted.pars = fit$par$`c(sd)` 
fitted.sd = round(fitted.pars[1], 3)
fitted.sd

fitted.vals = tibble(x=rnorm(n=10000, mean=0, sd=fitted.sd))
p2 <- p + 
  geom_density(data = fitted.vals, aes(x=x), color = 'blue', linetype = 'dashed')
p2
target_dir = here("data", "figs")
ggsave(paste(target_dir, "deviation_stochastic_ind.png", sep=FS), p2, 
       width = 8, height = 4)


# fit normal(mu,sd) hdi -------------------------------------------------
# highest density interval:
remove_extremes = 0.2
cutoff1 <- quantile(draws.ind %>% pull(diff), probs = remove_extremes/2)
cutoff2 <- quantile(draws.ind %>% pull(diff), probs = 1 - remove_extremes/2)

p + geom_vline(aes(xintercept = cutoff1)) + geom_vline(aes(xintercept = cutoff2))
ind.hdi = draws.ind %>% filter(diff >= cutoff1 & diff <= cutoff2)

mu <- greta::uniform(-0.1, 0.1, 1)
sd <- greta::uniform(0, 0.1, 1) 
greta::distribution(ind.hdi$diff) <- greta::normal(mu, sd)
model = greta::model(c(mu, sd))

fit = opt(model)
fitted.pars = fit$par$`c(mu, sd)` 
fitted.mu = round(fitted.pars[1], 3)
fitted.sd = round(fitted.pars[2], 3)

fitted.mu
fitted.sd
fitted.vals = tibble(x=rnorm(n=1000, mean=fitted.mu, sd=fitted.sd))

p2 <- p + 
  geom_vline(aes(xintercept = cutoff1)) + geom_vline(aes(xintercept = cutoff2)) +
  geom_density(data = fitted.vals, aes(x=x), color = 'blue', linetype = 'dashed')
p2
fit

target_dir = here("data", "figs")
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = T)
ggsave(paste(target_dir, "deviation_stochastic_ind.png", sep=FS), p2, 
       width = 8, height = 4)


