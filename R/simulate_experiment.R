library(ggplot2)
library(dplyr)

simulate_experiment <- function(params, treatments, seed = NULL,
                                write_grid = FALSE, ramp_stress = FALSE,
                                out_dir = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    set.seed(sample.int(1e9, 1))
  }

  # derived params
  params$effects <- list(
    drought = params$effect_frac$drought * params$mu,
    heat = params$effect_frac$heat * params$mu,
    drought_heat = params$effect_frac$drought_heat * params$mu
  )

  # Generate all plants and store in a dataframe where each row is a unique
  # combination of run x chamber x treatment x plant
  plants <- expand.grid(
    run = seq_len(params$n_runs),
    chamber = seq_len(params$n_chambers),
    treat = treatments,
    plant = seq_len(params$n_plants),
    KEEP.OUT.ATTRS = FALSE # prevents metadata related to levels, as messes with mutate
  ) 
  
  # Make a readable plant ID and add to the dataframe, make treatment a factor
  plants <- plants %>%
    mutate(
      treat = factor(treat, levels = treatments),
      plant_id = sprintf("r%d_c%d_%s_p%d", run, chamber, treat, plant)
    )

  # Sub-sample measured plants - N plants per run × chamber × treatment
  measured_plants <- plants %>%
    group_by(run, chamber, treat) %>%
    slice_sample(n = params$n_measured) %>%
    ungroup()

  # Check number of sampled plants
  # should be 3 runs * 2 chambers * 4 treatments * 6 plants = 144
  stopifnot(
    nrow(measured_plants) ==
      params$n_runs * params$n_chambers * length(treatments) * params$n_measured
  )

  # Create the full design, keep note of which plants we're measuring
  plants_full <- plants %>%
    arrange(run, chamber, treat, plant) %>%
    group_by(run, chamber, treat) %>%
    mutate(tray = ceiling(plant / params$n_measured)) %>%
    ungroup() %>%
    mutate(measured = plant_id %in% measured_plants$plant_id) %>%
    select(run, chamber, treat, plant, tray, plant_id, measured)

  if (write_grid) {
    write.csv(plants_full,
              file = file.path(out_dir, "full_experiment_grid.csv"),
              row.names = FALSE)
  }
  
  # Create plant x week combinations
  plant_week <- expand.grid(
    plant_id = measured_plants$plant_id,
    week = seq_len(params$n_weeks),
    KEEP.OUT.ATTRS = FALSE
  )
  
  # Add in metadata to the overal experimental design grid
  design <- plant_week %>%
    left_join(measured_plants, by = "plant_id")
  
  # Add factors for modelling
  design <- design %>%
    mutate(
      treat = factor(treat, levels = treatments),
      drought = as.integer(treat %in% c("drought", "heat_drought")),
      temp    = as.integer(treat %in% c("heat", "heat_drought"))
    )
  

  # simulate random effects

  # run effects
  rand_eff_run <- rnorm(params$n_runs, 0, params$sd$run)
  names(rand_eff_run) <- as.character(seq_len(params$n_runs))

  # chamber within run effect
  run_chamber_levels <- with(design, interaction(run, chamber, drop = TRUE))
  rand_eff_chamber <- rnorm(length(unique(run_chamber_levels)), 0, params$sd$chamber)
  names(rand_eff_chamber) <- unique(run_chamber_levels)

  # plant effect
  rand_eff_plant <- rnorm(length(unique(design$plant_id)), 0, params$sd$plant)
  names(rand_eff_plant) <- unique(design$plant_id)



  # fixed effects

  stress_weeks <- (params$n_weeks - params$n_stress_weeks + 1):params$n_weeks

  effect_map <- tibble(
    treat = factor(treatments, levels = treatments),
    base_effect = c(
      control = 0,
      drought = -params$effects$drought,
      heat = -params$effects$heat,
      heat_drought = -params$effects$drought_heat
    )[treatments]
  )

  design <- design %>%
    left_join(effect_map, by = "treat")

  if (!ramp_stress) {
    design <- design %>%
      mutate(mu = params$mu + ifelse(week %in% stress_weeks, base_effect, 0))
  } else {
    design <- design %>%
      mutate(
        stress_progress = pmax(0, week - min(stress_weeks) + 1) /
          length(stress_weeks),
        mu = params$mu + base_effect * pmin(stress_progress, 1)
      )
  }

  # combine random + fixed effects
  design <- design %>%
    mutate(
      # run-to-run variability
      # Each run has a slightly different "baseline" Anet
      rand_eff_run = rand_eff_run[as.character(run)],

      # chamber-within-run variability
      # Each chamber within a run may differ slightly
      # Nested: run:chamber
      rand_eff_chamber = rand_eff_chamber[as.character(interaction(run, chamber))],

      # plant-to-plant variability
      # individual plants differ biologically even under same treatment
      rand_eff_plant = rand_eff_plant[as.character(plant_id)],

      # residual week-to-week variation
      # captures natural fluctuations in Anet over repeated measurements
      resid = rnorm(n(), 0, params$sd$resid),

      # combine fixed and random effects
      Anet = mu + rand_eff_run + rand_eff_chamber +
        rand_eff_plant + resid
    )

  list(
    design = design,
    plants_full = plants_full
  )
}
