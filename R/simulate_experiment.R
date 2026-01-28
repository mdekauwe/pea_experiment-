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

  # Convert fractional treatment effects into absolute effects
  params$effects <- list(
    drought = params$mu * params$effect_frac$drought,
    heat = params$mu * params$effect_frac$heat,
    drought_heat = params$mu * params$effect_frac$drought_heat
  )

  
  # Construct the experimental setup "table", with one row per unique
  # combination of run x chamber x plant
  plants_df <- expand.grid(
    run = seq_len(params$n_runs),
    chamber = seq_len(params$n_chambers),
    plant = seq_len(params$n_plants),
    KEEP.OUT.ATTRS = FALSE
  )
  
  # Define treatment assignment vectors for each chamber.
  # Each chamber receives a repeating sequence of treatments, one  per plant in 
  # that chamber.
  # - Chamber 1 alternates between "control" and "drought"
  # - Chamber 2 alternates between "heat" and "heat_drought"
  treat_ch1 <- rep(c("control", "drought"), length.out = params$n_plants)
  treat_ch2 <- rep(c("heat", "heat_drought"), length.out = params$n_plants)
  
  # Assign treatments to plants based on their chamber.
  # For each run and chamber, add a 'treat' column:
  # - Plants in chamber 1 receive treatments from treat_ch1
  # - Plants in chamber 2 receive treatments from treat_ch2
  plants_df <- plants_df %>%
    group_by(run, chamber) %>%
    mutate(
      treat = ifelse(chamber == 1, treat_ch1, treat_ch2)
    ) %>%
    ungroup()
  
  # Ensure treat is a factor with the correct order of levels.
  # Levels are set based on the treatments vector.
  plants_df <- plants_df %>%
    mutate(treat = factor(treat, levels = treatments))
  
  # create a unique plant ID
  plants_df <- plants_df %>%
    mutate(
      plant_id = sprintf("r%d_c%d_%s_p%d", run, chamber, treat, plant)
    )
  
  # Assign tray numbers to plants within each experimental unit.
  # Experimental unit depends on treatment:
  # * Drought treatments are applied at the plant level
  # * Heat treatments are applied at the chamber level
  # Plants are first ordered by run, chamber, treatment, and plant ID
  # Each tray holds 6 plants. The ceiling(plant / 6) operation assigns:
  # plants 1–6 -> tray 1; plants 7–12 -> tray 2; etc
  plant_grid <- plants_df %>%
    arrange(run, chamber, treat, plant) %>%
    group_by(run, chamber, treat) %>%
    mutate(tray = ceiling(plant / 6)) %>%
    ungroup()
  
  # Randomise tray positions within each run × chamber
  plant_grid <- plant_grid %>%
    group_by(run, chamber) %>%
    mutate(
      tray_pos = setNames(sample(unique(tray)), unique(tray))[as.character(tray)]
    ) %>%
    ungroup()
  
  # Randomly select one plant per tray to be measured
  measured_plants <- plant_grid %>%
    group_by(run, chamber, tray) %>%  # only by tray
    slice_sample(n = 1) %>%           # one plant per tray
    ungroup()
  
  # Check number of sampled plants
  # calculate expected number based on treatments present per chamber
  stopifnot(all(
    measured_plants %>%
      group_by(run, chamber, tray) %>%
      summarise(n = n()) %>%
      pull(n) == 1
  ))

  # Create the full experimental design, keep note of which plants 
  # we're measuring
  plant_grid <- plant_grid %>%
    mutate(measured = plant_id %in% measured_plants$plant_id)

  if (write_grid) {
    write.csv(plant_grid,
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
  experiment_df <- plant_week %>%
    left_join(measured_plants, by = "plant_id")
  
  # Add factors for modelling
  experiment_df <- experiment_df %>%
    mutate(
      treat = factor(treat, levels = treatments),
      drought = as.integer(treat %in% c("drought", "heat_drought")),
      temp    = as.integer(treat %in% c("heat", "heat_drought"))
    )
  

  # simulate random effects
  
  # Create a string ID for run × chamber
  experiment_df <- experiment_df %>%
    mutate(run_chamber_id = paste0("r", run, "_c", chamber))
  
  # run effects, generate N random numbers from a noraml distribution
  rand_eff_run <- rnorm(params$n_runs, 0, params$sd$run)
  # name the vector by the run ID to match to the correct random effect
  names(rand_eff_run) <- as.character(seq_len(params$n_runs))

  # chamber within run effect
  rand_eff_chamber <- rnorm(length(unique(experiment_df$run_chamber_id)), 0, params$sd$chamber)
  names(rand_eff_chamber) <- unique(experiment_df$run_chamber_id)

  # plant effect
  rand_eff_plant <- rnorm(length(unique(experiment_df$plant_id)), 0, params$sd$plant)
  names(rand_eff_plant) <- unique(experiment_df$plant_id)


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
  experiment_df <- experiment_df %>%
    left_join(effect_map, by = "treat")
  
  # Determine when stress starts
  stress_start <- params$n_weeks - params$n_stress_weeks + 1
  stress_weeks <- stress_start:params$n_weeks
  
  if (gradual_stress) {
    
    experiment_df <- experiment_df %>%
      mutate(
        mu = case_when(
          week < stress_start ~ params$mu,
          week >= stress_start & week <= params$n_weeks ~
            params$mu + pmin(0, base_effect * ((week - stress_start + 1) / params$n_stress_weeks))
        )
      )
    
  } else {
    
    # sudden stress applied only during stress weeks
    experiment_df <- experiment_df %>%
      mutate(
        mu = ifelse(week %in% stress_weeks, params$mu + base_effect, params$mu)
      )
  }

  # combine random + fixed effects and calculate the simulated response for
  # each plant x week in the experimental dataset
  experiment_df <- experiment_df %>%
    mutate(
      rand_eff_run = ifelse(week == 1, 0, rand_eff_run[as.character(run)]),
      rand_eff_chamber = ifelse(week == 1, 0, rand_eff_chamber[run_chamber_id]),
      rand_eff_plant = ifelse(week == 1, 0, rand_eff_plant[as.character(plant_id)]),
      resid = ifelse(week == 1, 0, rnorm(n(), 0, params$sd$resid)),
      Anet = mu + rand_eff_run + rand_eff_chamber + rand_eff_plant + resid
      
    )

  return(list(
    experiment_df = experiment_df,
    plant_grid = plant_grid
  ))
}
