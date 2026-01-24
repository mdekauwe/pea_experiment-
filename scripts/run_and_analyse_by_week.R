# Fit each week indvidually


library(lme4)
library(lmerTest)
library(simr)
library(ggplot2)
library(dplyr)
library(broom.mixed)
library(emmeans)


rm(list = ls())

setwd("/Users/xj21307/Desktop/peas")  
source("R/simulate_experiment.R")



################################################################################

#
## Set your experiment parameters
#

params <- list(
  n_runs = 3,           # number of experimental repeats
  n_chambers = 2,       # number of growth chambers
  n_plants = 36,        # number of plants within a treatment
  n_measured = 6,       # number of measured plants
  n_weeks  = 8,         # number of weeks of the experiment
  n_stress_weeks = 2,   # number of weeks with stress

  mu = 10,              # baseline photosynthetic rate, umol m-2 s-1

  effect_frac = list(
    drought = 0.15,     # 15% reduction
    heat = 0.20,        # 20% reduction
    drought_heat = 0.25 # 25% reduction
  ),

  sd = list(
    run = 0.8,          # between run variability, umol m-2 s-1
    chamber  = 0.5,     # between chamber variability, umol m-2 s-1
    plant = 2.0,        # between plant variability, umol m-2 s-1
    resid = 1.5         # residual noise, umol m-2 s-1
  )
)

treatments <- c("control", "drought", "heat", "heat_drought")

out_dir <- "~/Desktop/"

################################################################################


s <- simulate_experiment(
  params = params,
  treatments = treatments,
  seed = 123,
  write_grid = TRUE,
  gradual_stress = TRUE,
  out_dir = out_dir
)

design <- s$design
plants <- s$plants_full

#############################################
## Fit split-plot mixed model
#############################################

design <- design %>%
  mutate(week_f = factor(week))  # new factor variable

# capture hierarchical structure of the split-plot experiment
m_split_week <- lmer(Anet ~ treat * week_f +
                  (1 | run/chamber) + # whole-plot variability for heat, nested
                  (1 | run:chamber:plant_id), # sub-plot variability for drought
                  data = design)


tidy(m_split_week, effects = "fixed") %>%
  filter(grepl("week_f", term)) %>%
  arrange(term) %>%
  print(n = 28)

# Analysis of Anet using a split-plot linear mixed model, treating week as a 
# factor, showed that drought significantly reduced Anet starting in week 8 
# (drought × week 8: estimate = -2.00, t = -3.94), while heat significantly 
# reduced Anet in weeks 7 and 8 (heat × week 7: estimate = -1.25, t = -2.46; 
# heat × week 8: estimate = -2.75, t = -5.43). The combination of drought and 
# heat also showed significant effects in weeks 7 and 8 
# (drought × heat × week 7: estimate = -1.06, t = -2.10; 
# drought × heat × week 8: estimate = -2.44, t = -4.82), \
# indicating that the joint impact of drought and heat on Anet was 
# strongest in the later weeks of the experiment.


em <- emmeans(m_split_week, ~ treat | week_f)
pdf("~/Desktop/emmeans_plot.pdf", width = 8, height = 6)
plot(em)
dev.off()