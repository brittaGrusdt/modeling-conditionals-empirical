library(here)
library(rwebppl)
library(tidyverse)
library(truncnorm)
library(greta)
library(ExpDataWrangling)
library(ModelUtils)


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

pe.ind = data.pe %>% filter(str_detect(id, "independent")) %>% 
  group_by(prolific_id, id) %>% 
  transmute(pa = AC + `A-C`, pc = `AC` + `-AC`, p_ac.ind = pa * pc, 
            diff = abs(AC - p_ac.ind))

p = pe.ind %>% ggplot() + geom_density(aes(x = diff))


# fit normal distribution
mu <- greta::uniform(-1, 0.1, 1)
sd <- greta::uniform(0, 0.1, 1) 
#greta::distribution(pe.ind$diff) <- greta::normal(mu, sd, truncation = c(0, 1))
greta::distribution(pe.ind$diff) <- greta::normal(mu, sd)
model = greta::model(c(mu, sd))

fit = opt(model)
fitted.pars = fit$par$`c(mu, sd)` 
fitted.mu = round(fitted.pars[1], 3)
fitted.sd = round(fitted.pars[2], 3)

#fitted.vals = tibble(x=rtruncnorm(100, a=0, b=1, mean=fitted.mu, sd=fitted.sd))
fitted.vals = tibble(x=rnorm(n=100, mean=fitted.mu, sd=fitted.sd))

p + geom_density(data = fitted.vals, aes(x=x), color = 'red')
fit

fitted.mu
fitted.sd



