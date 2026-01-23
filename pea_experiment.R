library(lme4)
library(simr)
library(ggplot2)
library(dplyr)

set.seed(123)

rm(list = ls())

#
## Parameters
#
n_chambers <- 2
n_treat <- 4
n_plants <- 36
n_weeks <- 8
n_runs <- 3
n_measured <- 6
n_stress_weeks <- 2   # number of weeks with stress


mu <- 10                   # photosynthetic rate, umol m-2 s-1

# define effect sizes
effect_drought <- 0.15 * mu       # 15% reduction
effect_heat <- 0.20 * mu          # 20% reduction
effect_drought_heat <- 0.25 * mu  # 25% reduction

sd_run <- 0.8        # between run variability, umol m-2 s-1
sd_chamber <- 0.5    # between chamber variability, umol m-2 s-1
sd_plant <- 2.0      # between plant variability, umol m-2 s-1
sd_resid <- 1.5      # residual noise, umol m-2 s-1

treatments <- c("control", "drought", "heat", "heat_drought")

#
## Generate all plants (fully crossed design)
#
plants <- expand.grid(
  run = 1:n_runs,
  chamber = 1:n_chambers,
  treat = treatments,
  plant = 1:n_plants
)

# Make a readable plant ID
plants$plant_id <- with(plants, 
                        paste0("r", run, "_c", chamber, "_", treat, "_p", plant))

#
## Sub-sample measured plants - 6 plants per run × chamber × treatment
#

# Split the plant table by run × chamber × treatment
plant_groups <- split(plants, list(plants$run, plants$chamber, plants$treat), 
                      drop = TRUE)

# Create an empty list to store sampled plants
sampled_groups <- vector("list", length(plant_groups))
names(sampled_groups) <- names(plant_groups)

# Randomly select n_measured plants per group
for (i in seq_along(plant_groups)) {
  df <- plant_groups[[i]]  
  selected_rows <- sample(nrow(df), n_measured)
  sampled_groups[[i]] <- df[selected_rows, ]
}

measured_plants <- do.call(rbind, sampled_groups)

# Check number of sampled plants
# should be 3 runs * 2 chambers * 4 treatments * 6 plants = 144
nrow(measured_plants)  

#
## Expand across weeks
#
design <- expand.grid(plant_id = measured_plants$plant_id, week = 1:n_weeks)
design <- merge(design, measured_plants, by = "plant_id")

# create treatment indicators
design$drought <- ifelse(design$treat %in% c("drought", "heat_drought"), 1, 0)
design$temp    <- ifelse(design$treat %in% c("heat", "heat_drought"), 1, 0)

#
## simulate random effects
#

# run effects
rand_eff_run <- rnorm(n_runs, 0, sd_run)
names(rand_eff_run) <- as.character(1:n_runs)

# Chamber within run effect
run_chamber_levels <- with(design, interaction(run, chamber))
rand_eff_chamber <- rnorm(length(unique(run_chamber_levels)), 0, sd_chamber)
names(rand_eff_chamber) <- unique(run_chamber_levels)

# Plant effect
rand_eff_plant <- rnorm(length(unique(design$plant_id)), 0, sd_plant)
names(rand_eff_plant) <- unique(design$plant_id)

#
## fixed effects
#

baseline <- mu  # 10 umol m-2 s-1

# drought effect only for drought alone
drought_effect <- ifelse(design$treat == "drought", -effect_drought, 0)

# heat effect only for heat alone
heat_effect <- ifelse(design$treat == "heat", -effect_heat, 0)

# interaction effect only for combined stress
interaction_effect <- ifelse(design$treat == "heat_drought", -effect_drought_heat, 0)

# treatment applies only in the last N weeks
stress_weeks <- (n_weeks - n_stress_weeks + 1):n_weeks

design$mu <- baseline + ifelse(design$week %in% stress_weeks,
                               drought_effect + heat_effect + interaction_effect, 
                               0)

#
## random effects
#

# run-to-run variability
# Each run has a slightly different "baseline" Anet
design$rand_eff_run <- rand_eff_run[as.character(design$run)]

# chamber-within-run variability
# Each chamber within a run may differ slightly
# Nested: run:chamber
design$rand_eff_chamber <- rand_eff_chamber[as.character(interaction(design$run, design$chamber))]

# plant-to-plant variability
# individual plants differ biologically even under same treatment
design$rand_eff_plant <- rand_eff_plant[as.character(design$plant_id)]

# residual week-to-week variation
# captures natural fluctuations in Anet over repeated measurements
design$resid <- rnorm(nrow(design), 0, sd_resid)

# combine fixed and random effects
design$Anet <- (design$mu + design$rand_eff_run + design$rand_eff_chamber + 
                design$rand_eff_plant + design$resid)

m <- lmer(Anet ~ drought * temp * factor(week) + 
            (1 | run) + (1 | run:chamber) + (1 | plant_id), 
          data = design)

summary(m)

# this simulates a baseline anet of 9.23
# drought = 0.19
# temp = 0.95
# week = 0.05
# drought:temp interactions = -0.28
# drought:week, effect of drought over changing weeks = -0.24
# temp:week, effect of heat over changing weeks= -0.37
# drought:temp:week, three way interaction = 0.24


#
## Plot the experiment
#
summary_df <- design %>%
  group_by(week, treat) %>%
  summarise(mean_Anet = mean(Anet), se_Anet = sd(Anet) / sqrt(n()), 
            .groups = "drop")

ggplot(summary_df, aes(x = week, y = mean_Anet, color = treat, group = treat)) +
  geom_line(size = 1) +                          
  geom_point(size = 2) +                         
  geom_ribbon(aes(ymin = mean_Anet - se_Anet, 
                  ymax = mean_Anet + se_Anet, 
                  fill = treat), alpha = 0.2, color = NA) +       
  labs(
    x = "Week",
    y = expression(paste("Anet (", mu, "mol ", m^-2, s^-1, ")")),
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

#
## Power analysis
#

# more than 80% is enough

powerSim(m, test = fixed("drought:temp", "t"), nsim = 200)
# Result? we would detect an effect of drought x heat 100% of the time, 


powerSim(m, test = fixed("temp:week", "t"), nsim = 200)
# Result? we would detect an effect of heat 100% of the time, 

powerSim(m, test = fixed("drought:temp:week", "t"), nsim = 200)
# Result? we would detect an effect of drought x heat 82% of the time, 


