library(lme4)
library(simr)
library(ggplot2)
library(dplyr)

set.seed(123)

rm(list = ls())

############################################# Mess with this
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
effect_heat <- 0.2 * mu          # 20% reduction
effect_drought_heat <- 0.25 * mu  # 25% reduction

sd_run <- 0.8        # between run variability, umol m-2 s-1
sd_chamber <- 0.5    # between chamber variability, umol m-2 s-1
sd_plant <- 2.0      # between plant variability, umol m-2 s-1
sd_resid <- 1.5      # residual noise, umol m-2 s-1

############################################################################## 

treatments <- c("control", "drought", "heat", "heat_drought")

#
## Generate all plants 
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
## Sub-sample measured plants - N plants per run × chamber × treatment
#

# Split the plant table by run × chamber × treatment
plant_groups <- split(plants, list(plants$run, plants$chamber, plants$treat), 
                      drop = TRUE)

# Create an empty list to store sampled plants
sampled_groups <- vector("list", length(plant_groups))
names(sampled_groups) <- names(plant_groups)

# Randomly select n_measured plants per group
for (i in seq_along(plant_groups)) { 
  df <- plant_groups[[i]]  # extract the i-th group as a dataframe
  selected_rows <- sample(nrow(df), n_measured) # pick random row without replacement
  sampled_groups[[i]] <- df[selected_rows, ] # store them
}

measured_plants <- do.call(rbind, sampled_groups)

# Check number of sampled plants
# should be 3 runs * 2 chambers * 4 treatments * 6 plants = 144
nrow(measured_plants)  

#
## Output the full design, keep note of which plants we're measuring
#
plants$measured <- ifelse(plants$plant_id %in% measured_plants$plant_id, TRUE, 
                          FALSE)

plants <- plants %>%
  arrange(run, chamber, treat, plant) %>% # just clean up the order
  group_by(run, chamber, treat) %>%
  mutate(tray = ceiling(plant / 6)) %>%  # create a new column: tray = 1-6 within each treatment
  ungroup() %>% # remove the grouping
  mutate(measured = plant_id %in% measured_plants$plant_id) %>% # add the measured column
  select(run, chamber, treat, plant, tray, plant_id, measured, everything()) # reorder the cols

write.csv(plants, "full_experiment_grid.csv", row.names = FALSE)


#
## Expand the experiment across weeks
#
design <- expand.grid(plant_id = measured_plants$plant_id, week = 1:n_weeks)
design <- merge(design, measured_plants, by = "plant_id")

# create treatment indicators
design$drought <- ifelse(design$treat %in% c("drought", "heat_drought"), 1, 0)
design$temp <- ifelse(design$treat %in% c("heat", "heat_drought"), 1, 0)

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

# drought effect only 
drought_effect <- ifelse(design$treat == "drought", -effect_drought, 0)

# heat effect only 
heat_effect <- ifelse(design$treat == "heat", -effect_heat, 0)

# interaction effect 
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


# week is a numeric, so we have a linear effect across week
#m <- lmer(Anet ~ drought * temp * week + (1 | run) + (1 | run:chamber) + 
#            (1 | plant_id), data=design)

# we could treat week as a categorical, each week would be independent 
#m <- lmer(Anet ~ drought * temp * factor(week) + (1 | run) + 
#            (1 | run:chamber) + (1 | plant_id), data = design)

#summary(m)


# capture hierarchical structure of the split-plot experiment
m_split <- lmer(Anet ~ drought * temp * week +
                  (1 | run/chamber) +         # whole-plot variability for heat, nested
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


## Plot the experiment
#
summary_df <- design %>%
  group_by(week, treat) %>%
  summarise(mean_Anet = mean(Anet), se_Anet = sd(Anet) / sqrt(n()), 
            .groups = "drop")

p <- ggplot(summary_df, aes(x = week, y = mean_Anet, color = treat, group = treat)) +
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
ggsave(filename = "Anet_experiment.png", plot = p, width = 8, height = 6, 
       dpi = 300)

print(p)
#
## Power analysis
#

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


## old below, which is wrong


# weeks as a numeric
# the model estimates one slope for week

#powerSim(m, test = fixed("drought:week", "t"), nsim = 200)
#Result? we would detect an effect of drought x heat 99% of the time, 

#powerSim(m, test = fixed("temp:week", "t"), nsim = 200)
#Result? we would detect an effect of heat 100% of the time, 

#powerSim(m, test = fixed("drought:temp:week", "t"), nsim = 200)
#Result? we would detect an effect of drought x heat 98% of the time, 



###
# weeks as a factor...
# the model estimates one slope for week

# then we'd have to do something like...
#powerSim(m, test = fixed("drought:week2", "t"), nsim = 200)
#powerSim(m, test = fixed("drought:week3", "t"), nsim = 200)

#powerSim(m, test = fixed("drought:week", "f"), nsim = 200)
#powerSim(m, test = fixed("temp:week", "f"), nsim = 200)
#powerSim(m, test = fixed("drought:temp:week", "f"), nsim = 200)
