library(ggplot2)
library(dplyr)

simulate_experiment <- function(params, treatments, seed = NULL,
                                write_grid = FALSE, gradual_stress = FALSE,
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
    plant = seq_len(params$n_plants),
    KEEP.OUT.ATTRS = FALSE
  )
  
  # Compute treatment vectors for each chamber
  treat_ch1 <- rep(c("control", "drought"), length.out = params$n_plants)
  treat_ch2 <- rep(c("heat", "heat_drought"), length.out = params$n_plants)
  
  # Assign treatments 
  plants <- plants %>%
    group_by(run, chamber) %>%
    mutate(
      treat = ifelse(chamber == 1, treat_ch1, treat_ch2)
    ) %>%
    ungroup()
  
  # Convert treat to a factor with the correct levels
  plants <- plants %>%
    mutate(treat = factor(treat, levels = treatments))
  
  # create a unique plant ID
  plants <- plants %>%
    mutate(
      plant_id = sprintf("r%d_c%d_%s_p%d", run, chamber, treat, plant)
    )

  # Sub-sample measured plants - N plants per run × chamber × treatment
  measured_plants <- plants %>%
    group_by(run, chamber, treat) %>%
    slice_sample(n = params$n_measured) %>%
    ungroup()
  
  # Check number of sampled plants
  # calculate expected number based on treatments present per chamber
  expected_plants <- measured_plants %>%
    group_by(run, chamber) %>%
    summarise(n_treat = n_distinct(treat), .groups = "drop") %>%
    summarise(total = sum(n_treat * params$n_measured)) %>%
    pull(total)
  
  stopifnot(nrow(measured_plants) == expected_plants)

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
  
  # Add in metadata to the overall experimental design grid
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
  
  # Create a string ID for run × chamber
  design <- design %>%
    mutate(run_chamber_id = paste0("r", run, "_c", chamber))
  
  # run effects, generate N random numbers from a noraml distribution
  rand_eff_run <- rnorm(params$n_runs, 0, params$sd$run)
  # name the vector by the run ID to match to the correct random effect
  names(rand_eff_run) <- as.character(seq_len(params$n_runs))

  # chamber within run effect
  rand_eff_chamber <- rnorm(length(unique(design$run_chamber_id)), 0, params$sd$chamber)
  names(rand_eff_chamber) <- unique(design$run_chamber_id)

  # plant effect
  rand_eff_plant <- rnorm(length(unique(design$plant_id)), 0, params$sd$plant)
  names(rand_eff_plant) <- unique(design$plant_id)


  # fixed effects
  
  # work out when the stress begins, for ease I've forced this to be at the
  # end of the experiment
  
  # determine when stress starts
  stress_start <- params$n_weeks - params$n_stress_weeks + 1
  
  # weeks under stress
  stress_weeks <- stress_start:params$n_weeks
  
  # create a vector treatment effects
  treatment_effects <- c(
    control = 0,
    drought = -params$effects$drought,
    heat = -params$effects$heat,
    heat_drought = -params$effects$drought_heat
  )
  
  # Create a small table (tibble) of the treatment and base effects so we can 
  # join to the main experimental design, so each plant in the experiment knows 
  # its base effect
  effect_map <- tibble(
    treat = factor(treatments, levels = treatments),
    base_effect = treatment_effects[treatments]
  )
  
  # attach the base_effects to each row in the design, so the simulation knows
  # which treatement effect to apply to each plant in each week
  design <- design %>%
    left_join(effect_map, by = "treat")

  if (!gradual_stress) {
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

  # combine random + fixed effects and calculate the simulated response for
  # each plant x week in the experimental dataset
  design <- design %>%
    mutate(
      
      # run-to-run variability
      # Each run has a slightly different "baseline" Anet
      # rand_eff_run is a vector of random numbers, one per experimental run
      rand_eff_run = rand_eff_run[as.character(run)], #index by run ID

      # chamber-within-run variability
      # Each chamber within a run may differ slightly
      # Nested: run:chamber
      # rand_eff_chamber is a vector of random numbers, one per chamber
      rand_eff_chamber = rand_eff_chamber[run_chamber_id],
      
      # plant-to-plant variability
      # individual plants differ biologically even under same treatment
      # rand_eff_plant is a vector of random numbers, one per plant
      rand_eff_plant = rand_eff_plant[as.character(plant_id)],

      # residual week-to-week variation - measurement error v
      # captures natural fluctuations in Anet over repeated measurements
      resid = rnorm(n(), 0, params$sd$resid),

      # combine fixed and random effects
      Anet = mu + rand_eff_run + rand_eff_chamber + rand_eff_plant + resid
    )

  return(list(
    design = design,
    plants_full = plants_full
  ))
}
