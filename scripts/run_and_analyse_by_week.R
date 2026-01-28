library(lme4)
library(lmerTest)
library(simr)
library(ggplot2)
library(dplyr)
library(broom.mixed)
library(emmeans)


rm(list = ls())

setwd("/Users/xj21307/src/R/pea_experiment")  
source("R/simulate_experiment.R")

################################################################################

#
## Experiment parameters
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

s <- simulate_experiment(params = params, treatments = treatments, seed = 124, 
                         write_grid = TRUE, gradual_stress = TRUE, 
                         out_dir = out_dir)

plant_grid <- s$plant_grid        # full plant layout (all plants)
experiment_df <- s$experiment_df  # measured plants across weeks

#############################################
## Fit split-plot mixed model
#############################################

# Week effects are interpreted as deviations relative to the final pre-stress 
# week
experiment_df <- experiment_df %>%
  mutate(
    week = factor(week, levels = 1:params$n_weeks),
    week = relevel(week, ref = as.character(params$n_weeks - params$n_stress_weeks))
  )

# estimate week-specific treatment deviations relative to the reference week
# tray position is fixed per pnat across weeks and treated as a random effect
m_split_week <- lmer(Anet ~ drought * temp * week + 
                       (1 | run/chamber) +            # chamber (whole-plot) random effect
                       (1 | run:chamber:tray_pos) +   # tray-level random effect
                       (1 | run:chamber:plant_id),    # plant-level random effect
                     data = experiment_df)
# plant-specific temporal trajectories
# (1 + week | run:chamber:plant_id)

fixed_eff <- tidy(m_split_week, effects = "fixed")

# Treatment × week effects
fixed_eff %>%
  filter(grepl("^week", term)) %>%
  arrange(term) %>%
  print(n = Inf)

# Three-way interaction
fixed_eff %>%
  filter(grepl("drought:temp:week", term)) %>%
  arrange(term) %>%
  print(n = Inf)


em <- emmeans(m_split_week, ~ drought * temp | week)

pdf("~/Desktop/emmeans_plot.pdf", width = 8, height = 6)
plot(em)
dev.off()