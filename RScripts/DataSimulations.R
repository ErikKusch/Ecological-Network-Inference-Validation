#' ####################################################################### #
#' PROJECT: [InfVal; Demographic Simulation of Interacting Species]
#' CONTENTS:
#'  - Simulation Calls across parallel cluster
#'  DEPENDENCIES:
#'  -
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
message("Running Simulation Framework")

pblapply(1:(n_runs + 1), FUN = function(ITER) {
  ITER <- ITER - 1
  if (ITER == 0) {
    FNAME <- file.path(Dir.Data, paste0(RunName, "_DEMO.RData"))
    t_inter <- 0.1 # when to record data
    n_spec <- 5
  } else {
    FNAME <- file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData"))
  }

  if (file.exists(FNAME)) {
    1 + 1
  } else {
    seed <- ITER * 1000 + 42


    # creating a network of interacting species
    Network_igraph <- Sim.Network(
      n_spec = n_spec, # how many species
      NetworkType = NetworkType, # if links are symmetric ("Association") or asymmetric ("Interaction")
      Sparcity = Sparcity, # how many links are exactly 0. Bound between 0 and 1
      MaxStrength = MaxStrength, # what absolute value is the maximum link strength
      seed = seed # seed for randomness
    )
    if (RunName == "NoSpaceBinaryInterac" | RunName == "BinaryInterac") {
      igraph::E(Network_igraph)$weight <- sign(igraph::E(Network_igraph)$weight)
    }

    # creating vector of carrying capacity for species
    CarryingK_vec <- Sim.CarryingK(
      n_spec = n_spec, # how many species
      k_range = k_range, # lower and upper limit of carrying capacities
      seed = seed # seed for randomness
    )

    # environmental preferences / niches of species
    Niches_vec <- Sim.Niche(
      n_spec = n_spec, # how many species
      Env_range = Env_range, # lower and upper bound of environmental preferences/niches
      seed = seed # seed for randomness
    )
    if (RunName == "NoSpaceBinaryInterac") {
      namesVec <- names(Niches_vec)
      Niches_vec <- rep(2, n_spec)
      names(Niches_vec) <- namesVec
    }

    # creating an initial dataframe of individuals of species
    Initialise_df <- Sim.Initialise(
      n_spec = n_spec, # how many species
      n_individuals = n_individuals, # how many individuals in data frame
      n_mode = n_mode, # whether to to n_individuals for each species or in total
      Env_range = Env_range, # lower and upper bound of environmental preferences/niches
      Trait_means = Niches_vec, # species-specific niches
      Trait_sd = Trait_sd, # standard deviation around species-specific niches for drawing individual-specific niches
      seed = seed # seed for randomness
    )
    if (RunName == "RareSpecies") {
      RareSpecies_vec <- sample(unique(Initialise_df$Species), size = round(0.3 * n_spec), replace = FALSE)
      for (i in RareSpecies_vec) {
        Initialise_df <- Initialise_df[-which(Initialise_df$Species == i)[1:round(0.75 * length(which(Initialise_df$Species == i)))], ]
      }
      CarryingK_vec[names(CarryingK_vec) %in% RareSpecies_vec] <- round(CarryingK_vec[names(CarryingK_vec) %in% RareSpecies_vec] * 0.25)
    }

    # creating environment for simulation
    if (RunName == "NoSpaceBinaryInterac") {
      x_gradient <- function(x) 1 # function relating x-location to environmental value
      y_gradient <- function(y) 1 # function relating y-location to environmental value
    }
    Env_mat <- Sim.Space(
      x_range = Env_range, # lower and upper bound of environmental coordinates in x dimension
      y_range = Env_range, # lower and upper bound of environmental coordinates in y dimension
      ncol = Env_col, nrow = Env_row, # number of cells in y and x dimension respectively
      x_gradient = x_gradient, # function relating x-location to environmental value
      y_gradient = y_gradient # function relating y-location to environmental value
    )

    if (ITER == 0) {
      save(
        list =
          c("Env_mat", "Network_igraph", "CarryingK_vec"),
        file = file.path(Dir.Data, paste0(RunName, "_DEMO_Environment.RData"))
      )
    }

    # actual simulation
    SimResult <- Sim.Compute(
      # Demographic parameters
      d0 = d0,
      b0 = b0,
      k_vec = CarryingK_vec,
      ID_df = Initialise_df,

      # Spatial parameters
      env.xy = Env_mat,
      env.sd = Env_sd,
      mig.sd = migration,
      mig.top = migration_top,
      mig.trunc = migration_trunc,

      # Interaction parameters
      interac.maxdis = Effect_Dis,
      interac.igraph = Network_igraph,
      interac.scale = 1,

      # Simulation parameters
      Sim.t.max = t_max,
      Sim.t.inter = t_inter,
      seed = seed,
      verbose = verbose, # whether to print progress in time as current time
      RunName = RunName,
      writeFile = FALSE
    )

    save(list = c("SimResult", "Network_igraph", "CarryingK_vec", "Niches_vec"), file = FNAME)
  }
})
