library(lme4)
library(lmerTest)
library(simr)
library(ggplot2)
library(dplyr)

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

# capture hierarchical structure of the split-plot experiment
m_split <- lmer(Anet ~ drought * temp * week +
                  (1 | run/chamber) + # whole-plot variability for heat, nested
                  (1 | run:chamber:plant_id), # sub-plot variability for drought
                  data = design)

summary(m_split)

# baseline = 9.2
# drought = 0.19
# temp = 0.95
# week = 0.05
# drought:temp = -0.28 # not significant, not enough power to detect
# drought:week = -0.24
# temp:week = -0.37
# drought:temp:week = 0.24773 # significant,
# so total slope = 0.05678-0.24342-0.37146+ 0.24773 = -0.3
# week + drought:week + temp:week +drought:temp:week
#

# Analysis of Anet using a split-plot linear mixed model showed that drought
# significantly influenced the decline in Anet
# (drought × week: estimate = -0.24, t = -4.13), and heat also has a significant
# effect (temp × week: estimate = -0.37, t = -6.31). A significant three-way
# interaction between drought, heat, and week
# (drought × temp × week: estimate = 0.25, t = 2.98), indicated that the
# combined impact of drought and heat on Anet over time differed from the
# sum of their individual effects.


##############################################
## Plot the experiment
##############################################

summary_df <- design %>%
  group_by(week, treat) %>%
  summarise(mean_Anet = mean(Anet), se_Anet = sd(Anet) / sqrt(n()),
            .groups = "drop")

p <- ggplot(summary_df, aes(x = week, y = mean_Anet, color = treat, 
            group = treat)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_Anet - se_Anet,
                  ymax = mean_Anet + se_Anet,
                  fill = treat), alpha = 0.2, color = NA) +
    labs(
      x = "Week",
      y = expression(Anet~"("*mu*" mol "*m^-2*" "*s^-1*")"),
      color = "Treatment",
      fill = "Treatment") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank())
ggsave(filename = file.path(out_dir, "Anet_experiment.png"),
       plot = p, width = 8, height = 6,
       dpi = 300)

print(p)


#############################################
## Power analysis
##############################################

# more than 80% is enough

# Number of experimental units per factor:
# - drought: number of plants per chamber/run
m_split <- extend(m_split, along = "plant_id", n = length(unique(design$plant_id)))
powerSim(m_split, test = fixed("drought:week", "t"), nsim = 200)

# - heat: number of chambers per run
chamber_means <- design %>%
  group_by(run, chamber, temp, week) %>%
  summarise(Anet = mean(Anet), .groups = "drop")

m_chamber <- lmer(Anet ~ temp * week + (1 | run), data = chamber_means)
powerSim(m_chamber, test = fixed("temp:week", "t"), nsim = 200)
