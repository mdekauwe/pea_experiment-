
library(lme4)
library(simr)
library(dplyr)

# Experimental parameters
params <- list(
  n_runs = 3,            # number of experimental repeats
  n_chambers = 2,        # number of growth chambers
  n_plants = 36,         # plants per treatment per chamber
  n_weeks = 8,           # total weeks of experiment
  n_stress_weeks = 2,    # weeks under stress
  mu = 10,               # baseline photosynthetic rate
  effect_frac = list(
    drought = 0.15,      # 15% reduction
    heat = 0.20,         # 20% reduction
    drought_heat = 0.25  # 25% reduction
  ),
  sd = list(
    run = 0.8,
    chamber = 0.5,
    plant = 2.0,
    resid = 1.5
  )
)

treatments <- c("control", "drought", "heat", "heat_drought")

# Construct minimal dataset
runs <- params$n_runs
chambers <- params$n_chambers
plants <- params$n_plants
weeks <- params$n_weeks

df <- expand.grid(
  run = factor(1:runs),
  chamber = factor(1:chambers),
  plant_id = factor(1:plants),
  week = 1:weeks
)

# Assign treatments 
df$treat <- factor(rep(treatments, length.out = nrow(df)))
df$drought <- as.integer(df$treat %in% c("drought","heat_drought"))
df$temp <- as.integer(df$treat %in% c("heat","heat_drought"))


# Gradual stress 
stress_start <- params$n_weeks - params$n_stress_weeks + 1

# Base treatment effects (negative numbers for reductions)
base_effect <- c(
  control = 0,
  drought = -params$effect_frac$drought * params$mu,
  heat = -params$effect_frac$heat * params$mu,
  heat_drought = -params$effect_frac$drought_heat * params$mu
)

df <- df %>%
  mutate(base_effect = base_effect[as.character(treat)]) %>%
  mutate(
    mu = case_when(
      week < stress_start ~ params$mu,
      week >= stress_start & week <= params$n_weeks ~
        params$mu + pmin(0, base_effect * ((week - stress_start + 1) / params$n_stress_weeks))
    )
  )


df$Anet <- df$mu

m <- lmer(Anet ~ drought * temp * week + (1 | run/chamber) + 
            (1 | run:chamber:plant_id), data = df, REML = FALSE)

# Set random effect SDs
VarCorr(m)$`run` <- params$sd$run^2
VarCorr(m)$`run:chamber` <- params$sd$chamber^2
VarCorr(m)$`run:chamber:plant_id` <- params$sd$plant^2
sigma(m) <- params$sd$resid

summary(m)

#
## Power analysis
#

# Extend along plant_id to test larger sample sizes
m_ext <- extend(m, along = "plant_id", n = 50)

# Power for drought × week
power_drought <- powerSim(m_ext, test = fixed("drought:week", "t"), nsim = 200)
print(power_drought)

# Power for heat × week
m_ext_heat <- extend(m, along = "chamber", n = 6)  # 6 chambers per run, for example
power_heat <- powerSim(m_ext_heat, test = fixed("temp:week", "t"), nsim = 200)
print(power_heat)

# Power for three-way interaction
m_ext_threeway <- extend(m, along = "chamber", n = 6)  # e.g., 6 chambers
power_threeway <- powerSim(m_ext_threeway, test = fixed("drought:temp:week", "t"), nsim = 200)
print(power_threeway)
